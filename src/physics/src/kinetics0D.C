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
      _n_species(_species_vars.n_species()),
      _mu_index(0),
      _k_index(0),
      _cp_index(0),
      _p0(10000)
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

    // Gas Evaluator
    Evaluator gas_evaluator( *(this->_gas_mixture) );

    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Number> &F_T = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    // Will be 1 since we are dealing with a 0 dimensional problem in space
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // Temperature
    std::vector<libMesh::Real> T;
    T.resize(n_qpoints);

    // 2D vector even though n_qp = 1 to satisfy context.get_* needs
    std::vector< std::vector<libMesh::Real> > mass_fractions;
    mass_fractions.resize(n_qpoints);

    // Molecular weight of species (kg/kmol)
    std::vector<libMesh::Real> M;
    M.resize(n_qpoints);

    // Ideal Gas Constant
    std::vector<libMesh::Real> R;
    R.resize(n_qpoints);

    // Specific Heat
    std::vector<libMesh::Real> cp;
    cp.resize(n_qpoints);

    // Enthalpy
    std::vector<std::vector<libMesh::Real> > h_s;
    h_s.resize(n_qpoints);

    // Species Production Rates
    std::vector<std::vector<libMesh::Real> > omega_dot_s;
    omega_dot_s.resize(n_qpoints);

    // Convenience
    const VariableIndex s0_var = this->_species_vars.species(0);

    // The number of local degrees of freedom in each variable.
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi =
      context.get_element_fe(s0_var)->get_phi();

    // Assemble Residuals
    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {

        libMesh::Real T = context.interior_value(this->_temp_vars.T(), qp);
        libMesh::Real Tdot;
        context.interior_rate(this->_temp_vars.T(), qp, Tdot);

        libMesh::out <<"temp ... " << T  << std::endl;
        //libMesh::out <<"temp dot... " << Tdot  << std::endl;

        // Sanity check
        libmesh_assert_greater(T,0);

        mass_fractions[qp].resize(this->_n_species);
        h_s[qp].resize(this->_n_species);
        omega_dot_s[qp].resize(this->_n_species);

        // Various summations needed for residual evaluations
        libMesh::Real hwsum = 0;
        libMesh::Real xcsum = 0;
        libMesh::Real xsum = 0;
        libMesh::Real wdotsum = 0;

        // Pre-calculate above summations
        for( unsigned int s = 0; s < this->_n_species; s++ )
          {
            // Ensure species dont overshoot to negative values
            mass_fractions[qp][s] = std::max( 0.0, context.interior_value( this->_species_vars.species(s),qp ));

            h_s[qp][s] = gas_evaluator.h_s( T, s );
            cp[qp] = gas_evaluator.cp(T, _p0, mass_fractions[qp]);
            gas_evaluator.omega_dot( T, _p0, mass_fractions[qp], omega_dot_s[qp] );

            xsum += mass_fractions[qp][s];
            xcsum += mass_fractions[qp][s]*cp[qp];

            wdotsum += omega_dot_s[qp][s];
            hwsum += h_s[qp][s]*omega_dot_s[qp][s];

          }

        M[qp] = gas_evaluator.M_mix( mass_fractions[qp] );
        R[qp] = gas_evaluator.R_mix( mass_fractions[qp] );

        libMesh::out << "mass_frac: " <<mass_fractions[qp][0] << " "
                     << mass_fractions[qp][1] << " "  << mass_fractions[qp][2] << " "  << std::endl;

        libMesh::out <<"omega dot " << omega_dot_s[qp][0] << " " <<  omega_dot_s[qp][1] << " " <<
          omega_dot_s[qp][2] << " " <<std::endl;

        // Temperature residual
        for (unsigned int i = 0; i != n_T_dofs; ++i)
          {
            F_T(i) += -(hwsum/xcsum) * T_phi[qp][0];

            libMesh::out<< "hws: " << hwsum << std::endl;
            libMesh::out<< "xcs: " << xcsum << std::endl;
            //libMesh::out<< "Tphi: " << T_phi[qp][0] << std::endl; == 1
            libMesh::out<< "F_T(i): " << F_T(i) << std::endl;
          }

        // Species residuals
        for(unsigned int s=0; s < this->_n_species; s++)
          {
            libMesh::DenseSubVector<libMesh::Number> &F_s =
              context.get_elem_residual(this->_species_vars.species(s));

            for (unsigned int i = 0; i != n_s_dofs; ++i)
              {
                F_s(i) +=  omega_dot_s[qp][s]*s_phi[i][qp]
                  - (mass_fractions[qp][s] * (wdotsum/xsum +  Tdot/T ) )*s_phi[i][qp];
                //libMesh::out<< "F_s(i): " << F_s(i) << std::endl;
                //libMesh::out<< "s_phi: " << s_phi[i][qp] << std::endl; == 1
              }
          }
      }

if( compute_jacobian )
      libmesh_not_implemented();
  } // end element_time_derivative


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
                                     PhysicsNaming::kinetics0D(),
                                     (*this),
                                     _p0 );
  }

  template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::set_time_evolving_vars( libMesh::FEMSystem * system )
  {
    // Tell the system to march temperature forward in time
    system->time_evolving(_temp_vars.T(), 1);

    // And also each of the involved species
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
        system->time_evolving( _species_vars.species(i), 1 );
      }
}


} // namespace GRINS
