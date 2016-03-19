//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_FACTORY_WITH_GETPOT_H
#define GRINS_FACTORY_WITH_GETPOT_H

// libMesh
#include "libmesh/factory.h"
#include "libmesh/getpot.h"

namespace GRINS
{
  //! Abstract factory that provides availability of GetPot
  template<typename Base>
  class FactoryWithGetPot : public libMesh::Factory<Base>
  {
  public:
    FactoryWithGetPot( const std::string& name )
      : libMesh::Factory<Base>(name)
    {}

    ~FactoryWithGetPot(){};

    static void set_getpot( const GetPot& input )
    { _input = &input; }

  protected:

    /*! We store only a raw pointer here because we *can't* make a copy.
        Otherwise, the UFO detection will be all screwed. We are not taking
        ownership of this, so we need to *not* delete this.*/
    static const GetPot* _input;

  };

} // end namespace GRINS

#endif // GRINS_FACTORY_WITH_GETPOT_H
