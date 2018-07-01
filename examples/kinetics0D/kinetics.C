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

// GRINS
#include "grins/runner.h"
#include "grins/string_utils.h"
#include "grins/simulation.h"
#include "grins/multiphysics_sys.h"
#include "grins/composite_qoi.h"

// libMesh
#include "libmesh/time_solver.h"
#include "libmesh/numeric_vector.h"


#if GRINS_HAVE_ANTIOCH
#include "grins/ignition_delay_qoi.h"
#include "grins/antioch_chemistry.h"
#endif

int main(int argc, char* argv[])
{
#if GRINS_HAVE_ANTIOCH

  // Initialize runner for the simulation
  GRINS::Runner runner(argc,argv);
  runner.init();

  // Parse the input file
  // const GetPot & input = runner.get_input_file();

  // Get some objects well need (for convenience)
  GRINS::Simulation & sim = runner.get_simulation();
  GRINS::MultiphysicsSystem * system = sim.get_multiphysics_system();

  // Run the simulation
  runner.run();

  // Now build a data structure of the solution to house their history
  //system->get_mesh().read("hydrogen.xda");
  //libMesh::out << "mesh read in!...." << std::endl;

  // Assemble and print the QoI info
  GRINS::CompositeQoI * comp_qoi = libMesh::cast_ptr<GRINS::CompositeQoI*>(system->get_qoi());
  // GRINS::IgnitionDelayQoI & ig_qoi = libMesh::cast_ref<GRINS::IgnitionDelayQoI &>(comp_qoi->get_qoi(0));

  libMesh::QoISet qs;
  qs.add_index(0);
  //system->assemble_qoi(qs);
  libMesh::Real qoi_value = sim.get_qoi_value(0);

  libMesh::out << "printing my qoi...." << std::endl;
  if (system->get_mesh().comm().rank() == 0)
    std::cout <<std::fixed <<std::setprecision(16) << qoi_value <<std::endl;

#else
  libmesh_error_msg("ERROR: GRINS must be built with Antioch to use this Kinetics0D example. Please reconfigure your build to include the Antioch library.");
#endif
  return 0;
}
