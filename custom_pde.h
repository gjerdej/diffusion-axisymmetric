// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/type_enums.h>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>

#include <random>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, dim, ScalarValue>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , D(get_user_inputs().user_constants.get_double("D"))
    , c0(get_user_inputs().user_constants.get_double("c0"))
    , A(get_user_inputs().user_constants.get_double("A"))
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    const dealii::Tensor<1, dim> &mesh_size = get_user_inputs().spatial_discretization.rectangular_mesh.size;

    if (index == 0)
      {
        const auto &mesh_size = get_user_inputs().spatial_discretization.rectangular_mesh.size;

        const double R = mesh_size[0];
        const double H = mesh_size[1];
        const double r = point[0];
        const double z = point[1];

        const double pi    = std::acos(-1.0);
        const double alpha = 3.831705970207512;

        scalar_value = c0 + A * boost::math::cyl_bessel_j(0.0, alpha * r / R) * std::cos(pi * z / H);
      }
  }

  void
  compute_rhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    const double dt = sim_timer.get_timestep();
    if (solve_block_id == 0) // c
      {
        const auto point = variable_list.get_q_point_location();
        const ScalarValue r = point[0];

        ScalarValue c_old = variable_list.template get_value<Scalar, OldOne>(0);

        ScalarValue eq_c = r * c_old;

        variable_list.set_value_term(0, eq_c);
      }
    else if (solve_block_id == 1) // analytical solution
      {
        const auto point = variable_list.get_q_point_location();

        const auto &mesh_size = get_user_inputs().spatial_discretization.rectangular_mesh.size;

        const double R = mesh_size[0];
        const double H = mesh_size[1];
        const double t = sim_timer.get_time();

        const double pi    = std::acos(-1.0);
        const double alpha = 3.831705970207512;

        const double k_squared = (alpha / R) * (alpha / R) + (pi / H) * (pi / H);

        ScalarValue c_exact;

        // Evaluate scalar special functions separately in each SIMD lane.
        for (unsigned int lane = 0; lane < ScalarValue::size(); ++lane)
          {
            const double r = point[0][lane];
            const double z = point[1][lane];

            c_exact[lane] = c0 + A * boost::math::cyl_bessel_j(0.0, alpha * r / R)
                              * std::cos(pi * z / H) * std::exp(-D[lane] * k_squared * t);
          }

        variable_list.set_value_term(1, c_exact);
      }
  }

  void
  compute_lhs([[maybe_unused]] FieldContainer<dim, degree, number> &variable_list,
              [[maybe_unused]] const SimulationTimer               &sim_timer,
              [[maybe_unused]] unsigned int solve_block_id) const override
  {
    const double dt = sim_timer.get_timestep();
    if (solve_block_id == 0) // c
      {
        const auto point = variable_list.get_q_point_location();
        const ScalarValue r = point[0];

        ScalarValue c   = variable_list.template get_value<Scalar, LHS>(0);
        ScalarGrad  cx  = variable_list.template get_gradient<Scalar, LHS>(0);

        ScalarValue eq_c = r * c;
        ScalarGrad eqx_c = r * sim_timer.get_timestep() * D * cx;

        variable_list.set_value_term(0, eq_c);
        variable_list.set_gradient_term(0, eqx_c);
      }
  }

  ScalarValue D;
  double c0;
  double A;
};

PRISMS_PF_END_NAMESPACE
