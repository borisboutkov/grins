//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_INCOMPRESSIBLE_HYPERELASTICITY_WEAK_FORM_H
#define GRINS_INCOMPRESSIBLE_HYPERELASTICITY_WEAK_FORM_H

#include "grins/hyperelasticity_weak_form.h"

namespace GRINS
{
  //! Class to encapsulate residual evaluations for hyperelasticity
  /*!
   * Several Physics classes will reuse these terms so let's encapsulate them
   * in an object. These will sit inside dof loops so we want them to be inlined.
   * This object supports both 2D (plane strain) and 3D evaluations. The functions
   * are differentiated by the number of residual vectors passed to the function (2 or 3)
   * and similarly for the Jacobian.
   */
  template<typename StrainEnergy>
  class IncompressibleHyperelasticityWeakForm : public HyperelasticityWeakForm<StrainEnergy>
  {
  public:

    IncompressibleHyperelasticityWeakForm() = default;
    virtual ~IncompressibleHyperelasticityWeakForm() = default;

  protected:

    //! Volumetric stored energy function
    /*!
     * We use the form U(J) = 1/2(ln(J)^2)
     */
    libMesh::Number U( libMesh::Number J ) const;

    //! First derivative of volumetric stored energy function
    libMesh::Number dU( libMesh::Number J ) const;

    //! Second derivative of volumetric stored energy function
    libMesh::Number d2U( libMesh::Number J ) const;

  };

  template<typename StrainEnergy>
  inline
  libMesh::Number IncompressibleHyperelasticityWeakForm<StrainEnergy>::U( libMesh::Number J ) const
  {
    libMesh::Number lnJ = std::log(J);
    return 0.5*lnJ*lnJ;
  }

  template<typename StrainEnergy>
  inline
  libMesh::Number IncompressibleHyperelasticityWeakForm<StrainEnergy>::dU( libMesh::Number J ) const
  {
    return std::log(J)/J;
  }

  template<typename StrainEnergy>
  inline
  libMesh::Number IncompressibleHyperelasticityWeakForm<StrainEnergy>::d2U( libMesh::Number J ) const
  {
    return (1.0-std::log(J))/(J*J);
  }

} // end namespace GRINS

#endif // GRINS_INCOMPRESSIBLE_HYPERELASTICITY_WEAK_FORM_H
