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
#include "grins/kinetics_0d.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/cached_quantities_enum.h"
#include "grins/generic_ic_handler.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  Kinetics0d::Kinetics0d( const std::string & physics_name, const GetPot & input )
    : ScalarODE( physics_name, input ),
      _species_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>(VariablesParsing::species_mass_frac_variable_name(input,physics_name,VariablesParsing::PHYSICS)))
  {
    this->_ic_handler = new GenericICHandler( physics_name, input );
  }


  void Kinetics0d::init_context( AssemblyContext& context )
  {

  }

  void Kinetics0d::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    if( compute_jacobian )
      libmesh_not_implemented();

    const CachedValues & cache = context.get_cached_values();

    // Convenience
    const VariableIndex s0_var = this->_species_vars.species(0);

    // The number of local degrees of freedom in each variable.
       const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();

    // Element Jacobian * quadrature weights for interior integration.
       const std::vector<libMesh::Real> JxW;

       // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi = context.get_element_fe(s0_var)->get_phi();

    // The species shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::Gradient> >& s_grad_phi = context.get_element_fe(s0_var)->get_dphi();


    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {

        const std::vector<libMesh::Real>& omega_dot =
          cache.get_cached_vector_values(Cache::OMEGA_DOT)[qp];

        const std::vector<libMesh::Real>& h =
          cache.get_cached_vector_values(Cache::SPECIES_ENTHALPY)[qp];

        libMesh::Real jac = JxW[qp];

        libMesh::Real M = cache.get_cached_values(Cache::MOLAR_MASS)[qp];

        std::vector<libMesh::Gradient> grad_ws = cache.get_cached_vector_gradient_values(Cache::MASS_FRACTIONS_GRAD)[qp];
        libmesh_assert_equal_to( grad_ws.size(), this->_n_species );

        // Species residuals
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            libMesh::DenseSubVector<libMesh::Number> &Fs =
              context.get_elem_residual(this->_species_vars.species(s)); // R_{s}

            const libMesh::Real term1 = -omega_dot[s];

            for (unsigned int i=0; i != n_s_dofs; i++)
              {
                /*! \todo Need to add SCEBD term. */
                Fs(i) += ( term1*s_phi[i][qp] )*jac;
              }
          }

      } // quadrature loop

    return;
  }

  void Kinetics0d::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    const VariableIndex s0_var = this->_species_vars.species(0);
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();

    const std::vector<std::vector<libMesh::Real> >& s_phi =
      context.get_element_fe(s0_var)->get_phi();

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    // Element Jacobian * quadrature weights for interior integration.
       const std::vector<libMesh::Real> JxW;


    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        std::vector<libMesh::Real> ws(this->_n_species);
        for(unsigned int s=0; s < this->_n_species; s++ )
          ws[s] = context.interior_value(this->_species_vars.species(s), qp);

        libMesh::Real jac = JxW[qp];

        // M_dot = -M^2 \sum_s w_dot[s]/Ms
        libMesh::Real M_dot = 0.0;

        // Species residual
        for(unsigned int s=0; s < this->_n_species; s++)
          {
            libMesh::DenseSubVector<libMesh::Number> &F_s =
              context.get_elem_residual(this->_species_vars.species(s));

            libMesh::Real ws_dot;
            context.interior_rate(this->_species_vars.species(s), qp, ws_dot);

            for (unsigned int i = 0; i != n_s_dofs; ++i)
              F_s(i) -= ws_dot*s_phi[i][qp]*jac;

            // Start accumulating M_dot
            M_dot += ws_dot;
          }



        if( compute_jacobian )
          libmesh_not_implemented();

      }
  }

} // namespace GRINS
