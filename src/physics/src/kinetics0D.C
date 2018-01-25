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

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{

  template<typename Mixture, typename Evaluator>
  Kinetics0D<Mixture,Evaluator>::Kinetics0D(const PhysicsName& physics_name,
                                            const GetPot& input,
                                            std::unique_ptr<Mixture> & gas_mix)
    : /*Kinetics0dBase<Mixture>(physics_name,input,gas_mix),*/
    _rho_index(0),
    _mu_index(0),
    _k_index(0),
    _cp_index(0)
  {
    this->_ic_handler = new GenericICHandler( physics_name, input );
  }

  template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::init_context( AssemblyContext& context )
  {
    // First call base class
    //Kinetics0dAbstract::init_context(context);

    // We also need the side shape functions, etc.
    context.get_side_fe(this->_temp_vars.T())->get_JxW();
    context.get_side_fe(this->_temp_vars.T())->get_phi();
    context.get_side_fe(this->_temp_vars.T())->get_dphi();
    context.get_side_fe(this->_temp_vars.T())->get_xyz();
  }


   template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    if( compute_jacobian )
      libmesh_not_implemented();

    const CachedValues & cache = context.get_cached_values();

    // Convenience
    const VariableIndex s0_var = this->_species_vars.species(0);

    // The number of local degrees of freedom in each variable.
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();
        const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();

        // Element Jacobian * quadrature weights for interior integration.
        const std::vector<libMesh::Real>& JxW;

        // The species shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi = context.get_element_fe(s0_var)->get_phi();

    // The species shape function gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::Gradient> >& s_grad_phi = context.get_element_fe(s0_var)->get_dphi();

    // The temperature shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    // The temperature shape functions gradients at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& T_gradphi =
      context.get_element_fe(this->_temp_vars.T())->get_dphi();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    libMesh::DenseSubVector<libMesh::Number> &FT = context.get_elem_residual(this->_temp_vars.T()); // R_{T}

    unsigned int n_qpoints = context.get_element_qrule().n_points();
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number rho = cache.get_cached_values(Cache::MIXTURE_DENSITY)[qp];

        libMesh::Number T;

        T = cache.get_cached_values(Cache::TEMPERATURE)[qp];

        const libMesh::Gradient& grad_T =
          cache.get_cached_gradient_values(Cache::TEMPERATURE_GRAD)[qp];


        libMesh::Number mu = cache.get_cached_values(Cache::MIXTURE_VISCOSITY)[qp];

        const std::vector<libMesh::Real>& D =
          cache.get_cached_vector_values(Cache::DIFFUSION_COEFFS)[qp];

        libMesh::Number cp =
          cache.get_cached_values(Cache::MIXTURE_SPECIFIC_HEAT_P)[qp];

        libMesh::Number k =
          cache.get_cached_values(Cache::MIXTURE_THERMAL_CONDUCTIVITY)[qp];

        const std::vector<libMesh::Real>& omega_dot =
          cache.get_cached_vector_values(Cache::OMEGA_DOT)[qp];

        const std::vector<libMesh::Real>& h =
          cache.get_cached_vector_values(Cache::SPECIES_ENTHALPY)[qp];


        libMesh::Real jac = JxW[qp];

        libMesh::Real M = cache.get_cached_values(Cache::MOLAR_MASS)[qp];

        std::vector<libMesh::Gradient> grad_ws = cache.get_cached_vector_gradient_values(Cache::MASS_FRACTIONS_GRAD)[qp];
        libmesh_assert_equal_to( grad_ws.size(), this->_n_species );

        // Continuity Residual
        libMesh::Gradient mass_term(0.0,0.0,0.0);
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            mass_term += grad_ws[s]/this->_gas_mixture->M(s);
          }
        mass_term *= M;

        // Species residuals
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            libMesh::DenseSubVector<libMesh::Number> &Fs =
              context.get_elem_residual(this->_species_vars.species(s)); // R_{s}

            const libMesh::Gradient term2 = -rho*D[s]*grad_ws[s];

            for (unsigned int i=0; i != n_s_dofs; i++)
              {
                /*! \todo Need to add SCEBD term. */
                Fs(i) += ( term2*s_grad_phi[i][qp] )*jac;
              }
          }

        // Energy residual
        libMesh::Real chem_term = 0.0;
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            chem_term += h[s]*omega_dot[s];
          }

        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            FT(i) += ( ( - chem_term )*T_phi[i][qp]
                       - k*grad_T*T_gradphi[i][qp]  )*jac;
          }

      } // quadrature loop

  }

  template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars.T()).size();
    const VariableIndex s0_var = this->_species_vars.species(0);
    const unsigned int n_s_dofs = context.get_dof_indices(s0_var).size();

    const std::vector<libMesh::Real> &JxW;

    const std::vector<std::vector<libMesh::Real> >& p_phi =
      context.get_element_fe(this->_press_var.p())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& T_phi =
      context.get_element_fe(this->_temp_vars.T())->get_phi();

    const std::vector<std::vector<libMesh::Real> >& s_phi =
      context.get_element_fe(s0_var)->get_phi();



    // The subvectors and submatrices we need to fill:
    libMesh::DenseSubVector<libMesh::Real> &F_T = context.get_elem_residual(this->_temp_vars.T());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    const std::vector<libMesh::Point>& u_qpoint =
      context.get_element_fe(this->_flow_vars.u())->get_xyz();

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {
        libMesh::Real T_dot;
        context.interior_rate(this->_temp_vars.T(), qp, T_dot);

        libMesh::Real T = context.interior_value(this->_temp_vars.T(), qp);

        std::vector<libMesh::Real> ws(this->n_species());
        for(unsigned int s=0; s < this->_n_species; s++ )
          ws[s] = context.interior_value(this->_species_vars.species(s), qp);

        Evaluator gas_evaluator( *(this->_gas_mixture) );
        const libMesh::Real R_mix = gas_evaluator.R_mix(ws);
        const libMesh::Real p0 = this->get_p0_steady(context,qp);
        const libMesh::Real rho = this->rho(T, p0, R_mix);
        const libMesh::Real cp = gas_evaluator.cp(T,p0,ws);
        const libMesh::Real M = gas_evaluator.M_mix( ws );

        libMesh::Real jac = JxW[qp];
        const libMesh::Number r = u_qpoint[qp](0);

        // M_dot = -M^2 \sum_s w_dot[s]/Ms
        libMesh::Real M_dot = 0.0;

        // Species residual
        for(unsigned int s=0; s < this->n_species(); s++)
          {
            libMesh::DenseSubVector<libMesh::Number> &F_s =
              context.get_elem_residual(this->_species_vars.species(s));

            libMesh::Real ws_dot;
            context.interior_rate(this->_species_vars.species(s), qp, ws_dot);

            for (unsigned int i = 0; i != n_s_dofs; ++i)
              F_s(i) -= rho*ws_dot*s_phi[i][qp]*jac;

            // Start accumulating M_dot
            M_dot += ws_dot/this->_gas_mixture->M(s);
          }

        // Continuity residual
        // M_dot = -M^2 \sum_s w_dot[s]/Ms
        libMesh::Real M_dot_over_M = M_dot*(-M);

        // Energy residual
        for (unsigned int i = 0; i != n_T_dofs; ++i)
          F_T(i) -= rho*cp*T_dot*T_phi[i][qp]*jac;

        if( compute_jacobian )
          libmesh_not_implemented();

      }
  }

   template<typename Mixture, typename Evaluator>
  void Kinetics0D<Mixture,Evaluator>::compute_element_time_derivative_cache
  ( AssemblyContext & context )
  {
    CachedValues & cache = context.get_cached_values();

    Evaluator gas_evaluator( *(this->_gas_mixture) );

    const unsigned int n_qpoints = context.get_element_qrule().n_points();

    std::vector<libMesh::Real> T;

    T.resize(n_qpoints);

    std::vector<libMesh::Gradient> grad_T;

    grad_T.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > mass_fractions;
    std::vector<std::vector<libMesh::Gradient> > grad_mass_fractions;
    mass_fractions.resize(n_qpoints);
    grad_mass_fractions.resize(n_qpoints);

    std::vector<libMesh::Real> M;
    M.resize(n_qpoints);

    std::vector<libMesh::Real> R;
    R.resize(n_qpoints);

    std::vector<libMesh::Real> rho;
    rho.resize(n_qpoints);

    std::vector<libMesh::Real> cp;
    cp.resize(n_qpoints);

    std::vector<libMesh::Real> mu;
    mu.resize(n_qpoints);

    std::vector<libMesh::Real> k;
    k.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > h_s;
    h_s.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > D_s;
    D_s.resize(n_qpoints);

    std::vector<std::vector<libMesh::Real> > omega_dot_s;
    omega_dot_s.resize(n_qpoints);

    for (unsigned int qp = 0; qp != n_qpoints; ++qp)
      {

        T[qp] = context.interior_value(this->_temp_vars.T(), qp);
        grad_T[qp] = context.interior_gradient(this->_temp_vars.T(), qp);

        mass_fractions[qp].resize(this->_n_species);
        grad_mass_fractions[qp].resize(this->_n_species);
        h_s[qp].resize(this->_n_species);

        for( unsigned int s = 0; s < this->_n_species; s++ )
          {
            /*! \todo Need to figure out something smarter for controling species
              that go slightly negative. */
            mass_fractions[qp][s] = std::max( context.interior_value(this->_species_vars.species(s),qp), 0.0 );
            grad_mass_fractions[qp][s] = context.interior_gradient(this->_species_vars.species(s),qp);
            h_s[qp][s] = gas_evaluator.h_s( T[qp], s );
          }

        M[qp] = gas_evaluator.M_mix( mass_fractions[qp] );

        R[qp] = gas_evaluator.R_mix( mass_fractions[qp] );

        rho[qp] = this->rho( T[qp], 0,/*p0[qp],*/ R[qp] );

        cp[qp] = gas_evaluator.cp(T[qp], 0/*p0[qp]*/, mass_fractions[qp]);

        D_s[qp].resize(this->_n_species);

        gas_evaluator.mu_and_k_and_D( T[qp], rho[qp], cp[qp], mass_fractions[qp],
                                      mu[qp], k[qp], D_s[qp] );

        omega_dot_s[qp].resize(this->_n_species);
        gas_evaluator.omega_dot( T[qp], rho[qp], mass_fractions[qp], omega_dot_s[qp] );
      }

    cache.set_values(Cache::TEMPERATURE, T);
    cache.set_gradient_values(Cache::TEMPERATURE_GRAD, grad_T);
    cache.set_vector_values(Cache::MASS_FRACTIONS, mass_fractions);
    cache.set_vector_gradient_values(Cache::MASS_FRACTIONS_GRAD, grad_mass_fractions);
    cache.set_values(Cache::MOLAR_MASS, M);
    cache.set_values(Cache::MIXTURE_GAS_CONSTANT, R);
    cache.set_values(Cache::MIXTURE_DENSITY, rho);
    cache.set_values(Cache::MIXTURE_SPECIFIC_HEAT_P, cp);
    cache.set_values(Cache::MIXTURE_VISCOSITY, mu);
    cache.set_values(Cache::MIXTURE_THERMAL_CONDUCTIVITY, k);
    cache.set_vector_values(Cache::DIFFUSION_COEFFS, D_s);
    cache.set_vector_values(Cache::SPECIES_ENTHALPY, h_s);
    cache.set_vector_values(Cache::OMEGA_DOT, omega_dot_s);
  }

} // namespace GRINS
