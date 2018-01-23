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


#ifndef GRINS_KINETICS_0D_H
#define GRINS_KINETICS_0D_H

#include "grins/scalar_ode.h"
#include "grins/multicomponent_variable.h"

namespace GRINS
{
  class Kinetics0d : public ScalarODE
  {
  public:

    Kinetics0d( const std::string & physics_name, const GetPot & input );

    virtual ~Kinetics0d(){}

    // Context initialization
    virtual void init_context( AssemblyContext& context );

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context );

    // Mass matrix part(s)
    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context );


  protected:

    //! Index from registering this quantity
    unsigned int _n_species;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _mole_fractions_index;

    //! Index from registering this quantity. Each species will have it's own index.
    std::vector<unsigned int> _omega_dot_index;

    SpeciesMassFractionsVariable & _species_vars;

  private:

    Kinetics0d();

  };

} // namespace GRINS

#endif //GRINS_KINETICS_0D_H
