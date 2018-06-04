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


#ifndef GRINS_KINETICS0D_H
#define GRINS_KINETICS0D_H

#include "grins/physics.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"

namespace GRINS
{
  template < typename Mixture, typename Evaluator >
    class Kinetics0D : public Physics  {

  public:

    // Constructor
    Kinetics0D(const PhysicsName& physics_name, const GetPot& input,
                                std::unique_ptr<Mixture> & gas_mix);

    // Destructor
    virtual ~Kinetics0D(){}

    // Context initialization
    virtual void init_context( AssemblyContext& context );

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context );

    // Sets temp variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Returns how many species are defined for this problem
    unsigned int n_species() const;

    // Function to evaluate gas mixture properties
    const Mixture & gas_mixture() const;

    // Evaluate rho
    libMesh::Real rho ( libMesh::Real T, libMesh::Real p0, libMesh::Real R_mix);

  protected:

    PrimitiveTempFEVariables & _temp_vars;

    SpeciesMassFractionsVariable & _species_vars;

    std::unique_ptr<Mixture> _gas_mixture;

    //! Index from registering this quantity
    unsigned int _n_species;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _species_viscosity;

    //! Index from registering this quantity
    unsigned int _mu_index;

    //! Index from registering this quantity
    unsigned int _k_index;

    //! Index from registering this quantity
    unsigned int _cp_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _mole_fractions_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _h_s_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _omega_dot_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _Ds_index;

    // Thermodynamic pressure (assumed constant)
    libMesh::Number _p0;


  private:

    // Read options from GetPot input file.
    void read_input_options( const GetPot& input );

    Kinetics0D();

  }; // class Kinetics0D

} // namespace GRINS

#endif //GRINS_KINETICS0D_H
