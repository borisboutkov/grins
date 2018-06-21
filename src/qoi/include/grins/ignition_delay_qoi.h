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


#ifndef GRINS_IGNITION_DELAY_QOI_H
#define GRINS_IGNITION_DELAY_QOI_H

// GRINS
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/auto_ptr.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

namespace GRINS
{
  //! Ignition Delay QoI
  /*!
    This class implements a QoI that calculates the amount of
    time until ignition of chemical mixture.
  */
  class IgnitionDelayQoI : public QoIBase
  {
  public:

    //! Constructor
    /*! Constructor takes GetPot object to read any input options associated
      with this QoI as well as a chemistry pointer */
    IgnitionDelayQoI( const GetPot & input, const std::string& qoi_name, const std::shared_ptr<GRINS::AntiochChemistry> chem );

    virtual ~IgnitionDelayQoI();

    //! Required to provide clone (deep-copy) for adding QoI object to libMesh objects.
    virtual QoIBase* clone() const;


    //! Override the QoIBase implementation
    virtual void parallel_op( const libMesh::Parallel::Communicator & communicator,
                              libMesh::Number & sys_qoi,
                              libMesh::Number & local_qoi );


    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    //! Initialize solution history info
    virtual void init( const GetPot& input,
                       const MultiphysicsSystem& system,
                       unsigned int qoi_num );


    //! Initialize context
    virtual void init_context( AssemblyContext& context );

    //! Compute the qoi value.
    virtual void element_qoi( AssemblyContext& context,
                              const unsigned int qoi_index );

    //! Compute the qoi derivative with respect to the solution.
    virtual void element_qoi_derivative( AssemblyContext& context,
                                         const unsigned int qoi_index );

  protected:

    //! Antioch/Cantera object
    std::shared_ptr<GRINS::AntiochChemistry> _chemistry;

    PrimitiveTempFEVariables & _T_var;
    SpeciesMassFractionsVariable & _Y_var;

    //! Reference temperature [K]
    libMesh::Real _T0;

    //! Index for the species of interest
    unsigned int _species_idx;

    //! Manual copy constructor due to the UniquePtr
    IgnitionDelayQoI(const IgnitionDelayQoI& original);

    //! Read in solution history from output files. Save to internal struct
    void read_solution_history( const GetPot& input, AssemblyContext& context );

    //! Storage for solution history post processing
    std::vector< std::vector<libMesh::Real> > _solution_history;

    //! Where to read solution data from
    std::string _filename_prefix;

    //! Which species is the fuel
    std::string _fuel_species;

    //! Percetage of initial fuel consumed
    std::string _fuel_consumption;

    //! Number of timesteps in the simulation
    unsigned int _n_timesteps;

  private:
    //! User never call default constructor.
    IgnitionDelayQoI();

    GetPot _input;

  };

  inline
  bool IgnitionDelayQoI::assemble_on_interior() const
  {
    return true;
  }

  inline
  bool IgnitionDelayQoI::assemble_on_sides() const
  {
    return false;
  }
}
#endif //GRINS_IGNITION_DELAY_QOI_H
