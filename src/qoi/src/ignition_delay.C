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
#include "grins/string_utils.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_base.h"
#include "libmesh/parsed_fem_function.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif


namespace GRINS
{
  IgnitionDelayQoI::IgnitionDelayQoI( const GetPot & input, const std::string& qoi_name, const std::shared_ptr<GRINS::AntiochChemistry>& chem )
    : QoIBase(qoi_name),
      _chemistry(chem),
      _T_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>("Temperature")),
      _Y_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>("SpeciesMassFractions")),
      _T0(296), // [K]
      _input(input)
  {

    // Parse input file for general IgnitionDelay parameters
    if( input.have_variable("vis-options/vis_output_file_prefix") )
      _filename_prefix = input("vis-options/vis_output_file_prefix","None");
    else
      libmesh_error_msg("ERROR: Could not QoI/IgnitionDelay/filename_prefix!");

    if( input.have_variable("SolverOptions/TimeStepping/n_timesteps") )
        _n_timesteps = input("SolverOptions/TimeStepping/n_timesteps",0);
    else
      libmesh_error_msg("ERROR: Could not find valid entry for n_timesteps!");

    if( input.have_variable("SolverOptions/TimeStepping/delta_t") )
        _dt = input("SolverOptions/TimeStepping/delta_t",0.0);
    else
      libmesh_error_msg("ERROR: Could not find valid entry for delta_t");

    // Ignition delay type specific input
    if( input.have_variable("QoI/IgnitionDelay/ignition_type") )
      {
        _ignition_type = input("QoI/IgnitionDelay/ignition_type","None");

        if ( _ignition_type == "fuel_consumption" )
          {
            if ( input.have_variable("QoI/IgnitionDelay/fuel_species") )
              _fuel_species = input("QoI/IgnitionDelay/fuel_species", "None");
            else
              libmesh_error_msg("ERROR: Could not QoI/IgnitionDelay/fuel_species!");

            if ( input.have_variable("QoI/IgnitionDelay/fuel_consumption") )
              _fuel_consumption = input("QoI/IgnitionDelay/fuel_consumption", 0.1);
            else
              libmesh_error_msg("ERROR: Could not find QoI/IgnitionDelay/fuel_consumption");
          }

        else if ( _ignition_type == "species_concentration_rate" )
          {
            if ( input.have_variable("QoI/IgnitionDelay/fuel_species") )
              _fuel_species = input("QoI/IgnitionDelay/fuel_species", "None");
            else
              libmesh_error_msg("ERROR: Could not find QoI/IgnitionDelay/fuel_species!");
          }
      }
    else
      libmesh_error_msg("ERROR: Could not find QoI/IgnitionDelay/ignition_type! Choose from : fuel_consumption, species_concentration_rate");



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


  void IgnitionDelayQoI::read_solution_history( const GetPot& input, AssemblyContext& context )
  {
    // Do the actual construction of the solution history
    libMesh::out << "reading solution history data..."<< std::endl;
    // libMesh cant read back in a 0-dim xda file... system calls for now!
    /*
      libMesh::Mesh temp_mesh (system.get_mesh().comm());
      temp_mesh.read(filename_prefix+".xda");
      libMesh::EquationSystems eqsys (temp_mesh);
      eqsys.init();
    */

    // size the _solution_history
    int n_species = context.get_dof_indices().size();
    this->_solution_history.resize(_n_timesteps, std::vector<libMesh::Real>(n_species));

    // Clean output data using the Kinetics0D/clean_kinetics.sh script... (See note above)
    std::string clean_data_command = "./clean_kinetics.sh " + _filename_prefix + " " + _filename_prefix+ ".in";
    std::system(clean_data_command.c_str());

    // Read the generated csv into _solution_history
    std::ifstream data("./" + _filename_prefix + "/clean_data.csv");

    int row = 0;
    std::string line;
    // Read away header info line and line of spaces
    std::getline(data,line);
    std::getline(data,line);

    // now parse the rest of the file
    while (std::getline(data,line))
      {
        std::stringstream  lineStream(line);
        std::vector<libMesh::Real> one_row;
        GRINS::StringUtilities::split_string_real( lineStream.str(), "," , one_row );
        _solution_history[row] = one_row;
        row++;
      }

    // check we have what we want
    /*
    libMesh::out << "printing soln_history data..."<< std::endl;
    libMesh::out << "of size: "<<_solution_history.size() << "x" << _solution_history[0].size() << std::endl;

    for (int i = 0; i < _solution_history.size(); i++)
      {
        for (int j = 0; j < _solution_history[i].size(); j++)
          {
            libMesh::out << _solution_history[i][j] << " ";
          }
        libMesh::out << std::endl;
      }
    */
  }

  void IgnitionDelayQoI::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* element_fe;
    context.get_element_fe<libMesh::Real>(0, element_fe);

  }

  void IgnitionDelayQoI::element_qoi( AssemblyContext& context,
                                       const unsigned int qoi_index )
  {
    this->read_solution_history(_input, context);

    /*! \todo Need to generalize this to the multiple QoI case */
    libMesh::Number& qoi = context.get_qois()[qoi_index];

    // + 1 because we added timestep data to the clean_data file
    unsigned int fuel_var_num = context.get_system().variable_number( _fuel_species ) + 1 ;

    // to extract the timestep
    libMesh::Real ts = 0;
    libMesh::Real initial_fuel = _solution_history[0][fuel_var_num];

    if ( this->_ignition_type == "fuel_consumption" )
      {
        // igntion delay = specified percentage of fuel is consumed
        for (unsigned int i = 0 ; i < this->_solution_history.size(); ++i )
          {
            if (  this->_solution_history[i][fuel_var_num]  <
                  (1 -_fuel_consumption) * initial_fuel )
              {
                ts = this->_solution_history[i][0];
                qoi = ts*this->_dt;
                return;
              }
          }
      }
    else if ( this->_ignition_type == "species_concentration_rate" )
      {
        // ignition delay =  species concentration rate maximum
        // just do centered diff in time on the solution
        std::vector<libMesh::Real> concentration_rate;
        concentration_rate.resize(this->_solution_history.size() - 2);
        for (unsigned int i = 1 ; i < _solution_history.size() - 1; i++ )
          {
            concentration_rate[i-1] = std::abs((this->_solution_history[i][fuel_var_num]
                          - this->_solution_history[i-1][fuel_var_num]))  / (2*this->_dt);
          }
        // ts is max of the rate vec
        auto maxrate = std::max_element(concentration_rate.begin(), concentration_rate.end());
        auto dist = std::distance(concentration_rate.begin(), maxrate);
        ts = this->_solution_history[dist][0];
        qoi = (ts+_dt/2)*this->_dt;
      }
    else
      {
        libmesh_error_msg("ERROR: QoI/IgnitionDelay/ignition_type not understood!");
      }

  } //end element_qoi

} //namespace GRINS
