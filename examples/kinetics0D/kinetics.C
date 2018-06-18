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

#include "libmesh/time_solver.h"
#include "libmesh/memory_solution_history.h"
#include "libmesh/numeric_vector.h"


#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

int main(int argc, char* argv[])
{
#if GRINS_HAVE_ANTIOCH
  GRINS::Runner runner(argc,argv);
  runner.init();

  // Parse the input file
  const GetPot & input = runner.get_input_file();

  GRINS::Simulation & sim = runner.get_simulation();
  GRINS::MultiphysicsSystem * system = sim.get_multiphysics_system();

  libMesh::MemorySolutionHistory  hist (*system);
  system->get_time_solver().set_solution_history(hist);
  //hist.set_overwrite_previously_stored(true);

  runner.run();

  system->get_time_solver().retrieve_timestep();

  //std::map<std::string, std::unique_ptr<libMesh::NumericVector<libMesh::Number>> saved_vectors;
  //  for ( auto it = hist->stored_solutions.begin(); it = hist->stored_solutions.end(); it++)
  //  std::cout <<std::fixed <<std::setprecision(16) << (*it)->first <<std::endl;



  GRINS::CompositeQoI * comp_qoi = libMesh::cast_ptr<GRINS::CompositeQoI*>(system->get_qoi());
  libMesh::QoISet qs;
  qs.add_index(0);
  system->assemble_qoi(qs);
  libMesh::Real qoi = sim.get_qoi_value(0);

  if (system->get_mesh().comm().rank() == 0)
    std::cout <<std::fixed <<std::setprecision(16) <<qoi <<std::endl;


#else
  libmesh_error_msg("ERROR: GRINS must be built with Antioch to use the Kinetics0D example. Please reconfigure your build to include the Antioch library.");
#endif
  return 0;
}
