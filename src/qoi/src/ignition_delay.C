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
#include "grins/ignition_delay_qoi.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_base.h"
#include "libmesh/parsed_fem_function.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif


namespace GRINS
{
  IgnitionDelayQoI::IgnitionDelayQoI( const std::string& qoi_name, const std::shared_ptr<GRINS::AntiochChemistry> chem )
    : QoIBase(qoi_name),
      _chemistry(chem),
      _T_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>("Temperature")),
      _Y_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>("SpeciesMassFractions")),
      _T0(296) // [K]
   {
     // read output files!
   }

  IgnitionDelayQoI::IgnitionDelayQoI( const IgnitionDelayQoI& original )
    : QoIBase(original),
      _chemistry(original._chemistry),
      _T_var(original._T_var),
      _Y_var(original._Y_var),
      _T0(296) // [K]
  {
    // do what in the copy?
  }

  IgnitionDelayQoI::~IgnitionDelayQoI() {}

  QoIBase* IgnitionDelayQoI::clone() const
  {
    return new IgnitionDelayQoI( *this );
  }

 void IgnitionDelayQoI::parallel_op( const libMesh::Parallel::Communicator & communicator,
                                             libMesh::Number & sys_qoi,
                                             libMesh::Number & local_qoi )
  {
    communicator.max(local_qoi);
    sys_qoi = local_qoi;
    QoIBase::_qoi_value = sys_qoi;
  }


  void IgnitionDelayQoI::init
  (const GetPot& input,
   const MultiphysicsSystem& system,
   unsigned int /*qoi_num*/ )
  {
    libMesh::ParsedFEMFunction<libMesh::Number> *qf
      (new libMesh::ParsedFEMFunction<libMesh::Number>
       (system, ""));
    this->qoi_functional.reset(qf);

    this->set_parameter(*qf, input,
                        "QoI/ParsedInterior/qoi_functional", "DIE!");
  }

  void IgnitionDelayQoI::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* element_fe;
    context.get_element_fe<libMesh::Real>(0, element_fe);
    element_fe->get_JxW();
    element_fe->get_xyz();

    qoi_functional->init_context(context);
  }

  void IgnitionDelayQoI::element_qoi( AssemblyContext& context,
                                       const unsigned int qoi_index )
  {

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    /*! \todo Need to generalize this to the multiple QoI case */
    libMesh::Number& qoi = context.get_qois()[qoi_index];

    for( unsigned int qp = 0; qp != n_qpoints; qp++ )
      {
        qoi = context.time;
      }
  }

  void IgnitionDelayQoI::element_qoi_derivative( AssemblyContext& context,
                                                  const unsigned int qoi_index )
  {
    libMesh::DenseSubVector<libMesh::Number> & dQdT  = context.get_qoi_derivatives(qoi_index, _T_var.T());

    std::vector<libMesh::Point> qp(1);
    //qp[0] = qp_ref;

    std::unique_ptr< libMesh::FEBase > T_fe = libMesh::FEGenericBase<libMesh::Real>::build(context.get_elem().dim(), context.get_element_fe(_T_var.T())->get_fe_type() );
    const std::vector<std::vector<libMesh::Real> > & T_phi =T_fe->get_phi();
    T_fe->reinit(&(context.get_elem()),&qp);

    std::unique_ptr< libMesh::FEBase > Ys_fe = libMesh::FEGenericBase<libMesh::Real>::build(context.get_elem().dim(), context.get_element_fe(_Y_var.species(_species_idx))->get_fe_type() );
    const std::vector<std::vector<libMesh::Real> > & Ys_phi = Ys_fe->get_phi();
    Ys_fe->reinit(&(context.get_elem()),&qp);

    libMesh::Real T; // temperature
    std::vector<libMesh::Real> Y(_chemistry->n_species()); // mass fractions

    context.point_value(_T_var.T(), qp[0], T); // [K]

    // all mass fractions needed to get M_mix
    for (unsigned int s=0; s<_chemistry->n_species(); s++)
      context.point_value(_Y_var.species(s), qp[0], Y[s]);

    // temperature derivatives
    for (unsigned int j=0; j<dQdT.size(); j++)
      dQdT(j) += T_phi[j][0];

    // mass fraction deriv for all species
    for (unsigned int s=0; s<Y.size(); ++s)
      {
        libMesh::DenseSubVector<libMesh::Number> & dQdYi = context.get_qoi_derivatives(qoi_index,_Y_var.species(s));
        for (unsigned int j=0; j<dQdYi.size(); j++)
          dQdYi(j) += Ys_phi[j][0];
      }
  }

} //namespace GRINS
