// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>

using namespace prismspf;

int
main(int argc, char *argv[])
{
  // Initialize MPI
  prismspf::MPIInitFinalize mpi_init(argc, argv);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);
  std::string     parameters_filename = cli_options.get_parameters_filename();

  constexpr unsigned int dim    = 2; // TODO change to 3 (original app)
  constexpr unsigned int degree = 1; // TODO change to 1 (original app)

  std::vector<FieldAttributes> fields = {FieldAttributes("c"), FieldAttributes("c_exact")};

  SolveBlock c_block;
  c_block.id            = 0;
  c_block.solve_type    = Linear;
  c_block.solve_timing  = Initialized;
  c_block.field_indices = {0};
  c_block.dependencies_rhs = make_dependency_set(fields, {"old_1(c)"});
  c_block.dependencies_lhs = make_dependency_set(fields, {"lhs(c)", "grad(lhs(c))"});

  SolveBlock exact_block;
  exact_block.id               = 1;
  exact_block.solve_type       = Explicit;
  exact_block.solve_timing     = PostProcess;
  exact_block.field_indices    = {1};
  exact_block.dependencies_rhs = {};
  exact_block.dependencies_lhs = {};

  std::vector<SolveBlock> solve_blocks({c_block, exact_block});

  UserInputParameters<dim>       user_inputs(parameters_filename);
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(fields,
                                       solve_blocks,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}