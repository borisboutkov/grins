//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-


// This class
#include "grins/kinetics0D.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/cached_quantities_enum.h"
#include "grins/generic_ic_handler.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  template<typename Mixture, typename Evaluator>
  Kinetics0D <Mixture,Evaluator>::Kinetics0D(const PhysicsName& physics_name,
                                            const GetPot& input,
                                            std::unique_ptr<Mixture> & gas_mix)
    : Physics(physics_name,input),
      _temp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _species_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>(VariablesParsing::species_mass_frac_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _gas_mixture(gas_mix.release()),
      _rho_index(0),
      _mu_index(0),
      _k_index(0),
      _n_species(_species_vars.n_species()),
      _p0(0),
      _cp_index(0)
  {
    this->_ic_handler = new GenericICHandler( physics_name, input );
    this->read_input_options(input);
  }

  template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::init_context( AssemblyContext& context )
  {

    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_species_vars.species(0))->get_JxW();
    context.get_element_fe(_species_vars.species(0))->get_phi();
    context.get_element_fe(_species_vars.species(0))->get_dphi();

    // We also need the temp functions
    context.get_element_fe(this->_temp_vars.T())->get_JxW();
    context.get_element_fe(this->_temp_vars.T())->get_phi();
    context.get_element_fe(this->_temp_vars.T())->get_dphi();

  }


  template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {

    // Our gas evaluator
    Evaluator gas_evaluator( *(this->_gas_mixture) );

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Number> &F_T = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    // Will be 1 since we are dealing with a 0 dimensional problem in space
    unsigned int n_qpoints = context.get_element_qrule().n_points();
    libMesh::out << "num  quadpoints = " << n_qpoints << std::endl;
    libMesh::out << "num  species = " << _n_species << std::endl;

    // Temperature
    std::vector<libMesh::Real> T;
    T.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > mass_fractions;
    mass_fractions.resize(n_qpoints);

    // Molecular weight of species (kg/kmol)
    std::vector<libMesh::Real> M;
    M.resize(n_qpoints);

    // Ideal Gas Constant
    std::vector<libMesh::Real> R;
    R.resize(n_qpoints);

    // density
    std::vector<libMesh::Real> rho;
    rho.resize(n_qpoints);

    // specific heat
    std::vector<libMesh::Real> cp;
    cp.resize(n_qpoints);

    // enthalpy
    std::vector<std::vector<libMesh::Real> > h_s;
    h_s.resize(n_qpoints);

    // species production rates
    std::vector<std::vector<libMesh::Real> > omega_dot_s;
    omega_dot_s.resize(n_qpoints);

    // Convenience
    const VariableIndex s0_var = this->_species_vars.species(0);

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();

    libMesh::out << "num  tdofs = " << n_T_dofs << std::endl;
    libMesh::out << "num  sdofs = " << n_s_dofs << std::endl;

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi =
      context.get_element_fe(s0_var)->get_phi();

    libMesh::out << "assemble... " << std::endl;
    // Assemble Residuals
    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {

        libMesh::Real T = context.interior_value(this->_temp_vars.T(), qp);

        libMesh::out <<"temp... " << T  << std::endl;


        mass_fractions[qp].resize(this->_n_species);
        h_s[qp].resize(this->_n_species);
        omega_dot_s[qp].resize(this->_n_species);

        libMesh::out <<"omega dot_s size " <<omega_dot_s[qp].size()  << std::endl;

        gas_evaluator.omega_dot( T, rho[qp], mass_fractions[qp], omega_dot_s[qp] );



        libMesh::Real T_dot;
        context.interior_rate(this->_temp_vars.T(), qp, T_dot);

        libMesh::Real hwsum = 0;
        libMesh::Real xcsum = 0;

        for( unsigned int s = 0; s < this->_n_species; s++ )
          {
            // ensure species dont go slightly negative
            mass_fractions[qp][s] = std::max( 0.0, context.interior_value( this->_species_vars.species(s),qp ));
            h_s[qp][s] = gas_evaluator.h_s( T, s );
            cp[qp] = gas_evaluator.cp(T, 0/*p0[qp]*/, mass_fractions[qp]);

            xcsum += mass_fractions[qp][s]*cp[qp];
            hwsum += h_s[qp][s]*omega_dot_s[qp][s];
          }

        M[qp] = gas_evaluator.M_mix( mass_fractions[qp] );
        R[qp] = gas_evaluator.R_mix( mass_fractions[qp] );
        const libMesh::Real p0 = this->_p0;

        // Temperature residual
        for (unsigned int i = 0; i != n_T_dofs; ++i)
          F_T(i) -= (hwsum/xcsum) * T_phi[qp][0];

        // Species residual
        for(unsigned int s=0; s < this->n_species(); s++)
          {
            libMesh::DenseSubVector<libMesh::Number> &F_s =
              context.get_elem_residual(this->_species_vars.species(s)); //R_{s}

            // only
            for (unsigned int i = 0; i != n_s_dofs; ++i)
              F_s(i) -=omega_dot_s[i][qp]*s_phi[i][qp];

          }


      }

if( compute_jacobian )
      libmesh_not_implemented();

  } // end nonlocal_time_derivative


  template<typename Mixture, typename Evaluator>
  unsigned int Kinetics0D<Mixture,Evaluator>::n_species() const
  {
    return _n_species;
  }

  template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::read_input_options( const GetPot& input )
  {
    // Read thermodynamic pressure info
    MaterialsParsing::read_property( input,
                                     "ThermodynamicPressure",
                                     PhysicsNaming::reacting_low_mach_navier_stokes(),
                                     (*this),
                                     _p0 );

  }

  template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::set_time_evolving_vars( libMesh::FEMSystem * system )
  {
    // do something?
  }

  template<typename Mixture, typename Evaluator>
  libMesh::Real Kinetics0D<Mixture,Evaluator>::rho( libMesh::Real T, libMesh::Real p0, libMesh::Real R_mix)
  {
    return _p0;
  }

} // namespace GRINS
