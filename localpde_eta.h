#ifndef LOCALPDE_
#define LOCALPDE_

#include <interfaces/pdeinterface.h>

using namespace DOpE;
using namespace dealii;
//using namespace std;

#if DEAL_II_VERSION_GTE(9,3,0)
template<
  template<bool DH, typename VECTOR, int dealdim> class EDC,
  template<bool DH, typename VECTOR, int dealdim> class FDC,
  bool DH, typename VECTOR, int dealdim>
class LocalPDE : public PDEInterface<EDC, FDC, DH, VECTOR, dealdim>
#else
template<
  template<template<int, int> class DH, typename VECTOR, int dealdim> class EDC,
  template<template<int, int> class DH, typename VECTOR, int dealdim> class FDC,
  template<int, int> class DH, typename VECTOR, int dealdim>
class LocalPDE : public PDEInterface<EDC, FDC, DH, VECTOR, dealdim>
#endif
{
public:
  static void declare_params(ParameterReader &param_reader)
  {
    param_reader.SetSubsection("main parameters");
    param_reader.declare_entry("end_time", "0.0", Patterns::Double(0));
    param_reader.declare_entry("time intervals", "30", Patterns::Integer(0));

    param_reader.SetSubsection("Local PDE parameters");
    param_reader.declare_entry("constant_k", "0.0", Patterns::Double(0));
    param_reader.declare_entry("alpha_eps", "0.0", Patterns::Double(0));
    param_reader.declare_entry("gamma_penal", "0.0", Patterns::Double(0));
    param_reader.declare_entry("G_c", "0.0", Patterns::Double(0));
    param_reader.declare_entry("E_modulus", "0.0", Patterns::Double(0));
    param_reader.declare_entry("poisson_ratio_nu", "0.0", Patterns::Double(0));
    param_reader.declare_entry("eta", "0.0", Patterns::Double(0));

    param_reader.SetSubsection("Local Functional parameters");
    param_reader.declare_entry("alpha_tikhonov", "0.0", Patterns::Double(0));
    param_reader.SetSubsection("Force parameters");
    param_reader.declare_entry("force_increment", "0.0", Patterns::Double(0));
  }

  void update_gamma(double gamma_penal)
  {
    gamma_penal_ = gamma_penal;
  }

  void update_alpha_tik(double alpha_tik)
  {
    alpha_tik_ = alpha_tik;
  }

  LocalPDE(ParameterReader &param_reader) : state_block_component_(3, 0)
  {
    state_block_component_[2] = 1;

    param_reader.SetSubsection("main parameters");
    n_time_intervals_ = param_reader.get_integer("time_intervals");
    end_time_ = param_reader.get_double("end_time");

    param_reader.SetSubsection("Local PDE parameters");
    constant_k_ = param_reader.get_double("constant_k");
    alpha_eps_ = param_reader.get_double("alpha_eps");
    gamma_penal_ = param_reader.get_double("gamma_penal");
    G_c_ = param_reader.get_double("G_c");
    E_modulus_ = param_reader.get_double("E_modulus");
    poisson_ratio_nu_ = param_reader.get_double("poisson_ratio_nu");
    eta_ = param_reader.get_double("eta");

    param_reader.SetSubsection("Local Functional parameters");
    alpha_tik_ = param_reader.get_double("alpha_tikhonov");
    lame_coefficient_mu_ = E_modulus_ / (2.0 * (1.0 + poisson_ratio_nu_));
    lame_coefficient_lambda_ = (2.0 * poisson_ratio_nu_ * lame_coefficient_mu_) / (1.0 - 2 * poisson_ratio_nu_);

    param_reader.SetSubsection("Force parameters");
    force_increment_ = param_reader.get_double("force_increment");
  }

  void ElementEquation(const EDC<DH, VECTOR, dealdim> &edc, dealii::Vector<double> &local_vector, double scale, double)
  {
    assert(this->problem_type_ == "state");

    const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
    unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
    unsigned int n_q_points = edc.GetNQPoints();

    uvalues_.resize(n_q_points, Vector<double>(3));
    ugrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
    last_timestep_uvalues_.resize(n_q_points, Vector<double>(3));
    edc.GetValuesState("last_newton_solution", uvalues_);
    edc.GetGradsState("last_newton_solution", ugrads_);
    edc.GetValuesState("last_time_solution", last_timestep_uvalues_);

    const FEValuesExtractors::Vector displacement(0);
    const FEValuesExtractors::Scalar phase(2);

    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
    {
      Tensor<1, 2> u;
      Init_Displacement(u, uvalues_[q_point]);
      double pf = uvalues_[q_point](2);

      Tensor<2, 2> grad_u;
      Init_Displacement_Gradient(grad_u, ugrads_[q_point]);

      Tensor<1, 2> grad_pf;
      Init_Phasefield_Gradient(grad_pf, ugrads_[q_point]);

      double old_timestep_pf = last_timestep_uvalues_[q_point](2);

      double g = (1.0 - constant_k_) * pf * pf + constant_k_;

      double pen = 0.;
      if (pf > old_timestep_pf)
      {
        pen = gamma_penal_;
      }

      Tensor<2, 2> stress_term = grad_u;
      Apply_Elasticity_Tensor(stress_term);

      for (unsigned int i = 0; i < n_dofs_per_element; i++)
      {
        const Tensor<2, 2> phi_i_grads_u = state_fe_values[displacement].symmetric_gradient(i, q_point);
        const double phi_i_pf = state_fe_values[phase].value(i, q_point);
        const Tensor<1, 2> phi_i_grads_pf = state_fe_values[phase].gradient(i, q_point);

        local_vector(i) += scale * ((n_time_intervals_ / end_time_) * (pen * (pf - old_timestep_pf) 
        * phi_i_pf + eta_ * (pf - old_timestep_pf) * phi_i_pf)) * state_fe_values.JxW(q_point);
        local_vector(i) += scale * (g * scalar_product(stress_term, phi_i_grads_u) 
        + G_c_ * alpha_eps_ * grad_pf * phi_i_grads_pf - G_c_ / alpha_eps_ * (1.0 - pf) * phi_i_pf 
        + (1.0 - constant_k_) * scalar_product(stress_term, grad_u) * pf * phi_i_pf) * state_fe_values.JxW(q_point);
      }
    }
  }

  void ElementMatrix(const EDC<DH, VECTOR, dealdim> &edc, FullMatrix<double> &local_matrix, double scale, double)
  {
    const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
    unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
    unsigned int n_q_points = edc.GetNQPoints();

    uvalues_.resize(n_q_points, Vector<double>(3));
    ugrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));

    if (this->GetTime() > 0.)
    {
      if (this->problem_type_ == "state")
      {
        edc.GetValuesState("last_newton_solution", uvalues_);
        edc.GetGradsState("last_newton_solution", ugrads_);
        edc.GetValuesState("last_time_solution", last_timestep_uvalues_);
      }
      else
      {
        edc.GetValuesState("state", uvalues_);
        edc.GetGradsState("state", ugrads_);
        edc.GetValuesState("state_i-1", last_timestep_uvalues_);
        if (this->problem_type_ == "adjoint" || this->problem_type_ == "adjoint_hessian")
        {
          edc.GetValuesState("state_i+1", next_timestep_uvalues_);
        }
      }

      const FEValuesExtractors::Vector displacement(0);
      const FEValuesExtractors::Scalar phase(2);

      std::vector<Tensor<1, 2>> phi_u(n_dofs_per_element);
      std::vector<Tensor<2, 2>> phi_grads_u(n_dofs_per_element);
      std::vector<double> phi_pf(n_dofs_per_element);
      std::vector<Tensor<1, 2>> phi_grads_pf(n_dofs_per_element);

      Tensor<2, 2> Identity;
      Identity[0][0] = 1.0;
      Identity[1][1] = 1.0;

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int k = 0; k < n_dofs_per_element; k++)
        {
          phi_u[k] = state_fe_values[displacement].value(k, q_point);
          phi_grads_u[k] = state_fe_values[displacement].symmetric_gradient(k, q_point);
          phi_pf[k] = state_fe_values[phase].value(k, q_point);
          phi_grads_pf[k] = state_fe_values[phase].gradient(k, q_point);
        }

        Tensor<1, 2> u;
        Init_Displacement(u, uvalues_[q_point]);

        double pf = uvalues_[q_point](2);

        Tensor<2, 2> grad_u;
        Init_Displacement_Gradient(grad_u, ugrads_[q_point]);

        Tensor<1, 2> grad_pf;
        Init_Phasefield_Gradient(grad_pf, ugrads_[q_point]);

        double g = (1.0 - constant_k_) * pf * pf + constant_k_;
        double pen = 0.0;

        if (this->problem_type_ == "adjoint" || this->problem_type_ == "adjoint_hessian")
        {
          double pf_plus = next_timestep_uvalues_[q_point](2);
          if (pf_plus > pf)
          {
            pen = gamma_penal_;
          }
        }
        else
        {
          double old_timestep_pf = last_timestep_uvalues_[q_point](2);
          if (pf > old_timestep_pf)
          {
            pen = gamma_penal_;
          }
        }

        Tensor<2, 2> stress_term = grad_u;
        Apply_Elasticity_Tensor(stress_term);

        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          Tensor<2, 2> stress_term_LinU = phi_grads_u[i];
          Apply_Elasticity_Tensor(stress_term_LinU);

          for (unsigned int j = 0; j < n_dofs_per_element; j++)
          {
            local_matrix(i, j) += scale * ((n_time_intervals_ / end_time_) * (pen * phi_pf[i] * phi_pf[j] 
            + eta_ * phi_pf[i] * phi_pf[j]) + 2.0 * pf * (1.0 - constant_k_) * phi_pf[i] * scalar_product(stress_term, phi_grads_u[j]) 
            + g * scalar_product(stress_term_LinU, phi_grads_u[j])) * state_fe_values.JxW(q_point);

            local_matrix(i, j) += scale * (2.0 * pf * (1.0 - constant_k_) * phi_pf[j] * scalar_product(stress_term_LinU, grad_u) 
            + (1.0 - constant_k_) * scalar_product(stress_term, grad_u) * phi_pf[i] * phi_pf[j] 
            + G_c_ / alpha_eps_ * phi_pf[i] * phi_pf[j] + G_c_ * alpha_eps_ * phi_grads_pf[i] * phi_grads_pf[j]) * state_fe_values.JxW(q_point);
          }
        }
      }
    }
    else
    {
      assert(this->GetTime() == 0.);
      assert(this->problem_type_ == "adjoint" || this->problem_type_ == "adjoint_hessian");

      const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
      unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
      unsigned int n_q_points = edc.GetNQPoints();

      const FEValuesExtractors::Vector displacement(0);
      const FEValuesExtractors::Scalar phase(2);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          const Tensor<1, 2> phi_i_u = state_fe_values[displacement].value(i, q_point);
          const double phi_i_pf = state_fe_values[phase].value(i, q_point);
          for (unsigned int j = 0; j < n_dofs_per_element; j++)
          {
            const Tensor<1, 2> phi_j_u = state_fe_values[displacement].value(j, q_point);
            const double phi_j_pf = state_fe_values[phase].value(j, q_point);
            local_matrix(i, j) += scale * (phi_i_u * phi_j_u + phi_i_pf * phi_j_pf) * state_fe_values.JxW(q_point);
          }
        }
      }
    }
  }

  void ElementEquation_U(const EDC<DH, VECTOR, dealdim> &edc, dealii::Vector<double> &local_vector, double scale, double)
  {
    assert(this->problem_type_ == "adjoint");

    if (this->GetTime() == 0.)
    {
      const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
      unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
      unsigned int n_q_points = edc.GetNQPoints();

      zvalues_.resize(n_q_points, Vector<double>(3));
      last_timestep_zvalues_.resize(n_q_points, Vector<double>(3));
      edc.GetValuesState("last_newton_solution", zvalues_);
      edc.GetValuesState("last_time_solution", last_timestep_zvalues_);

      const FEValuesExtractors::Vector displacement(0);
      const FEValuesExtractors::Scalar phase(2);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        Tensor<1, 2> z_v;
        Init_Displacement(z_v, zvalues_[q_point]);

        double z_pf = zvalues_[q_point](2);

        Tensor<1, 2> last_timestep_z_v;
        Init_Displacement(z_v, last_timestep_zvalues_[q_point]);

        double last_timestep_z_pf = last_timestep_zvalues_[q_point](2);

        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          const Tensor<1, 2> phi_i_u = state_fe_values[displacement].value(i, q_point);
          const double phi_i_pf = state_fe_values[phase].value(i, q_point);

          local_vector(i) += scale * (phi_i_u * (z_v - last_timestep_z_v) 
          + phi_i_pf * (z_pf - last_timestep_z_pf)) * state_fe_values.JxW(q_point);
          ;
        }
      }
    }
    else
    {
      const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
      unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
      unsigned int n_q_points = edc.GetNQPoints();

      zvalues_.resize(n_q_points, Vector<double>(3));
      edc.GetValuesState("last_newton_solution", zvalues_);

      zgrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
      edc.GetGradsState("last_newton_solution", zgrads_);

      uvalues_.resize(n_q_points, Vector<double>(3));
      edc.GetValuesState("state", uvalues_);

      ugrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
      edc.GetGradsState("state", ugrads_);

      last_timestep_zvalues_.resize(n_q_points, Vector<double>(3));
      edc.GetValuesState("last_time_solution", last_timestep_zvalues_);

      if (this->GetTime() < end_time_)
      {
        next_timestep_uvalues_.resize(n_q_points, Vector<double>(3));
        edc.GetValuesState("state_i+1", next_timestep_uvalues_);
      }
      const FEValuesExtractors::Vector displacement(0);
      const FEValuesExtractors::Scalar phase(2);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        Tensor<1, 2> z_v;
        Init_Displacement(z_v, zvalues_[q_point]);

        Tensor<2, 2> z_vgrad;
        Init_Displacement_Gradient(z_vgrad, zgrads_[q_point]);

        double z_pf = zvalues_[q_point](2);

        Tensor<1, 2> z_pfgrad;
        Init_Phasefield_Gradient(z_pfgrad, zgrads_[q_point]);

        Tensor<1, 2> u;
        Init_Displacement(u, uvalues_[q_point]);

        Tensor<2, 2> grad_u;
        Init_Displacement_Gradient(grad_u, ugrads_[q_point]);

        double pf = uvalues_[q_point](2);
        double last_timestep_z_pf = last_timestep_zvalues_[q_point](2);

        double g = (1.0 - constant_k_) * pf * pf + constant_k_;
        double pen = 0.0;

        double t_M_term = z_pf;
        if (this->GetTime() < end_time_)
        {
          double pf_plus = next_timestep_uvalues_[q_point](2);
          if (pf_plus > pf)
          {
            pen = gamma_penal_;
          }
          t_M_term = 0.0;
        }

        Tensor<2, 2> stress_term = grad_u;
        Apply_Elasticity_Tensor(stress_term);

        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          const Tensor<2, 2> phi_grads_u = state_fe_values[displacement].symmetric_gradient(i, q_point);
          const double phi_pf = state_fe_values[phase].value(i, q_point);

          const Tensor<1, 2> phi_grads_pf = state_fe_values[phase].gradient(i, q_point);
          Tensor<2, 2> stress_term_LinU = phi_grads_u;
          Apply_Elasticity_Tensor(stress_term_LinU);

          local_vector(i) += scale * (t_M_term + (n_time_intervals_ / end_time_) * (pen * (z_pf - last_timestep_z_pf) * phi_pf 
          + eta_ * (z_pf - last_timestep_z_pf) * phi_pf) + g * scalar_product(stress_term_LinU, z_vgrad) 
          + 2.0 * pf * (1.0 - constant_k_) * phi_pf * scalar_product(stress_term, z_vgrad)) * state_fe_values.JxW(q_point);

          local_vector(i) += scale * (+G_c_ / alpha_eps_ * phi_pf * z_pf + G_c_ * alpha_eps_ * phi_grads_pf * z_pfgrad 
          + 2.0 * pf * (1.0 - constant_k_) * scalar_product(stress_term_LinU, grad_u) * z_pf 
          + (1 - constant_k_) * scalar_product(stress_term, grad_u) * phi_pf * z_pf) * state_fe_values.JxW(q_point);
        }
      }
    }
  }

  void ElementEquation_UT(const EDC<DH, VECTOR, dealdim> &edc, dealii::Vector<double> &local_vector, double scale, double)
  {
    assert(this->problem_type_ == "tangent");

    const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
    unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
    unsigned int n_q_points = edc.GetNQPoints();

    duvalues_.resize(n_q_points, Vector<double>(3));
    dugrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
    uvalues_.resize(n_q_points, Vector<double>(3));
    ugrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
    last_timestep_uvalues_.resize(n_q_points, Vector<double>(3));
    last_timestep_duvalues_.resize(n_q_points, Vector<double>(3));
    edc.GetValuesState("last_newton_solution", duvalues_);
    edc.GetGradsState("last_newton_solution", dugrads_);
    edc.GetValuesState("state", uvalues_);
    edc.GetGradsState("state", ugrads_);
    edc.GetValuesState("state_i-1", last_timestep_uvalues_);
    edc.GetValuesState("last_time_solution", last_timestep_duvalues_);

    const FEValuesExtractors::Vector displacement(0);
    const FEValuesExtractors::Scalar phase(2);

    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
    {
      Tensor<1, 2> du_v;
      Init_Displacement(du_v, duvalues_[q_point]);
     
      double du_pf = duvalues_[q_point](2);

      Tensor<2, 2> du_vgrad;
      Init_Displacement_Gradient(du_vgrad, dugrads_[q_point]);

      Tensor<1, 2> du_pfgrad;
      Init_Phasefield_Gradient(du_pfgrad, dugrads_[q_point]);

      Tensor<1, 2> u;
      Init_Displacement(u, uvalues_[q_point]);

      double pf = uvalues_[q_point](2);

      Tensor<2, 2> grad_u;
      Init_Displacement_Gradient(grad_u, ugrads_[q_point]);

      double du_pf_minus = last_timestep_duvalues_[q_point](2);
      double pf_minus = last_timestep_uvalues_[q_point](2);

      double g = (1.0 - constant_k_) * pf * pf + constant_k_;

      double pen = 0.0;
      if (pf > pf_minus)
      {
        pen = gamma_penal_;
      }

      Tensor<2, 2> stress_term = grad_u;
      Apply_Elasticity_Tensor(stress_term);

      for (unsigned int i = 0; i < n_dofs_per_element; i++)
      {
        const Tensor<2, 2> phi_grads_u = state_fe_values[displacement].symmetric_gradient(i, q_point);
        const double phi_pf = state_fe_values[phase].value(i, q_point);
        const Tensor<1, 2> phi_grads_pf = state_fe_values[phase].gradient(i, q_point);

        Tensor<2, 2> stress_term_LinU = phi_grads_u;
        ;
        Apply_Elasticity_Tensor(stress_term_LinU);
        Tensor<2, 2> stress_term_du = du_vgrad;
        Apply_Elasticity_Tensor(stress_term_du);

        local_vector(i) += scale * ((n_time_intervals_ / end_time_) * (pen * (du_pf - du_pf_minus) * phi_pf 
        + eta_ * (du_pf - du_pf_minus) * phi_pf) + g * scalar_product(stress_term_du, phi_grads_u) 
        + 2.0 * pf * (1.0 - constant_k_) * du_pf * scalar_product(stress_term, phi_grads_u)) * state_fe_values.JxW(q_point);

        local_vector(i) += scale * (2.0 * pf * (1.0 - constant_k_) * scalar_product(stress_term_du, grad_u) * phi_pf 
        + (1 - constant_k_) * du_pf * scalar_product(stress_term, grad_u) * phi_pf 
        + G_c_ / alpha_eps_ * phi_pf * du_pf + G_c_ * alpha_eps_ * phi_grads_pf * du_pfgrad) * state_fe_values.JxW(q_point);
      }
    }
  }

  void ElementEquation_UTT(const EDC<DH, VECTOR, dealdim> &edc, dealii::Vector<double> &local_vector, double scale, double /*scale_ico*/)
  {
    assert(this->problem_type_ == "adjoint_hessian");

    if (this->GetTime() == 0.)
    {
      const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
      unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
      unsigned int n_q_points = edc.GetNQPoints();

      dzvalues_.resize(n_q_points, Vector<double>(3));
      last_timestep_dzvalues_.resize(n_q_points, Vector<double>(3));
      edc.GetValuesState("last_newton_solution", dzvalues_);
      edc.GetValuesState("last_time_solution", last_timestep_dzvalues_);

      const FEValuesExtractors::Vector displacement(0);
      const FEValuesExtractors::Scalar phase(2);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        Tensor<1, 2> dz_v;
        Init_Displacement(dz_v, zvalues_[q_point]);

        double dz_pf = dzvalues_[q_point](2);

        Tensor<1, 2> last_timestep_dz_v;
        Init_Displacement(dz_v, last_timestep_zvalues_[q_point]);

        double last_timestep_dz_pf = last_timestep_dzvalues_[q_point](2);

        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          const Tensor<1, 2> phi_i_u = state_fe_values[displacement].value(i, q_point);
          const double phi_i_pf = state_fe_values[phase].value(i, q_point);
          local_vector(i) += scale * (phi_i_u * (dz_v - last_timestep_dz_v) 
          + phi_i_pf * (dz_pf - last_timestep_dz_pf)) * state_fe_values.JxW(q_point);
          ;
        }
      }
    }
    else
    {
      const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
      unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
      unsigned int n_q_points = edc.GetNQPoints();

      dzvalues_.resize(n_q_points, Vector<double>(3));
      dzgrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
      uvalues_.resize(n_q_points, Vector<double>(3));
      ugrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
      last_timestep_dzvalues_.resize(n_q_points, Vector<double>(3));
      edc.GetValuesState("last_newton_solution", dzvalues_);
      edc.GetGradsState("last_newton_solution", dzgrads_);
      edc.GetValuesState("state", uvalues_);
      edc.GetGradsState("state", ugrads_);
      edc.GetValuesState("last_time_solution", last_timestep_dzvalues_);

      if (this->GetTime() < end_time_)
      {
        next_timestep_uvalues_.resize(n_q_points, Vector<double>(3));
        edc.GetValuesState("state_i+1", next_timestep_uvalues_);
      }

      const FEValuesExtractors::Vector displacement(0);
      const FEValuesExtractors::Scalar phase(2);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        Tensor<1, 2> dz_v;
        Init_Displacement(dz_v, dzvalues_[q_point]);
        double dz_pf = dzvalues_[q_point](2);

        Tensor<2, 2> dz_vgrad;
        Init_Displacement_Gradient(dz_vgrad, dzgrads_[q_point]);

        Tensor<1, 2> dz_pfgrad;
        Init_Phasefield_Gradient(dz_pfgrad, dzgrads_[q_point]);

        Tensor<1, 2> u;
        Init_Displacement(u, uvalues_[q_point]);

        double pf = uvalues_[q_point](2);

        Tensor<2, 2> grad_u;
        Init_Displacement_Gradient(grad_u, ugrads_[q_point]);

        double last_timestep_dz_pf = dzvalues_[q_point](2);

        double g = (1.0 - constant_k_) * pf * pf + constant_k_;

        double pen = 0.0;
        double t_M_term = dz_pf; 
        if (this->GetTime() < end_time_)
        {
          double pf_plus = next_timestep_uvalues_[q_point](2);
          if (pf_plus > pf)
          {
            pen = gamma_penal_;
          }
          t_M_term = 0.0;
        }

        Tensor<2, 2> stress_term = grad_u;
        Apply_Elasticity_Tensor(stress_term);

        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          const Tensor<2, 2> phi_grads_u = state_fe_values[displacement].symmetric_gradient(i, q_point);
          const double phi_pf = state_fe_values[phase].value(i, q_point);
          const Tensor<1, 2> phi_grads_pf = state_fe_values[phase].gradient(i, q_point);

          Tensor<2, 2> stress_term_LinU = phi_grads_u;
          Apply_Elasticity_Tensor(stress_term_LinU);

          local_vector(i) += scale * (t_M_term + (n_time_intervals_ / end_time_) * (pen * (dz_pf - last_timestep_dz_pf) * phi_pf 
          + eta_ * (dz_pf - last_timestep_dz_pf) * phi_pf) + g * scalar_product(stress_term_LinU, dz_vgrad) 
          + (1.0 - constant_k_) * 2.0 * scalar_product(stress_term, dz_vgrad) * pf * phi_pf) * state_fe_values.JxW(q_point);

          local_vector(i) += scale * ((1.0 - constant_k_) * 2.0 * scalar_product(stress_term, phi_grads_u) * pf * dz_pf 
          + (1 - constant_k_) * scalar_product(stress_term, grad_u) * phi_pf * dz_pf + G_c_ / alpha_eps_ * phi_pf * dz_pf 
          + G_c_ * alpha_eps_ * phi_grads_pf * dz_pfgrad) * state_fe_values.JxW(q_point);
        }
      }
    }
  }


  void ElementEquation_UU(const EDC<DH, VECTOR, dealdim> &edc, dealii::Vector<double> &local_vector, double scale, double)
  {
    assert(this->problem_type_ == "adjoint_hessian");

    if (this->GetTime() == 0.)
    {
      unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
      unsigned int n_q_points = edc.GetNQPoints();

      const FEValuesExtractors::Vector displacement(0);
      const FEValuesExtractors::Scalar phase(2);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          local_vector(i) += 0.0;
        }
      }
    }
    else
    {
      const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
      unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
      unsigned int n_q_points = edc.GetNQPoints();

      uvalues_.resize(n_q_points, Vector<double>(3));
      ugrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
      zvalues_.resize(n_q_points, Vector<double>(3));
      zgrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
      duvalues_.resize(n_q_points, Vector<double>(3));
      dugrads_.resize(n_q_points, std::vector<Tensor<1, 2>>(3));
      edc.GetValuesState("state", uvalues_);
      edc.GetGradsState("state", ugrads_);
      edc.GetValuesState("last_newton_solution", zvalues_);
      edc.GetGradsState("last_newton_solution", zgrads_);
      edc.GetValuesState("last_newton_solution", duvalues_);
      edc.GetGradsState("last_newton_solution", dugrads_);

      const FEValuesExtractors::Vector displacement(0);
      const FEValuesExtractors::Scalar phase(2);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        Tensor<1, 2> u;
        Init_Displacement(u, uvalues_[q_point]);

        double pf = uvalues_[q_point](2);

        Tensor<2, 2> grad_u;
        Init_Displacement_Gradient(grad_u, ugrads_[q_point]);

        Tensor<1, 2> z_v;
        Init_Displacement(z_v, zvalues_[q_point]);

        double z_pf = zvalues_[q_point](2);

        Tensor<2, 2> z_vgrad;
        Init_Displacement_Gradient(z_vgrad, zgrads_[q_point]);

        Tensor<1, 2> z_pfgrad;
        Init_Phasefield_Gradient(z_pfgrad, zgrads_[q_point]);

        Tensor<1, 2> du_v;
        Init_Displacement(du_v, duvalues_[q_point]);

        Tensor<2, 2> du_vgrad;
        Init_Displacement_Gradient(du_vgrad, dugrads_[q_point]);

        double du_pf = duvalues_[q_point](2);

        Tensor<1, 2> du_pfgrad;
        Init_Phasefield_Gradient(du_pfgrad, dugrads_[q_point]);

        Tensor<2, 2> stress_term = grad_u;
        Apply_Elasticity_Tensor(stress_term);

        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          const Tensor<2, 2> phi_grads_u = state_fe_values[displacement].symmetric_gradient(i, q_point);

          const double phi_pf = state_fe_values[phase].value(i, q_point);

          Tensor<2, 2> stress_term_LinU = phi_grads_u;
          Apply_Elasticity_Tensor(stress_term_LinU);

          Tensor<2, 2> stress_term_du = du_vgrad;
          Apply_Elasticity_Tensor(stress_term_du);


          local_vector(i) += 0.0 * scale * (2.0 * pf * (1.0 - constant_k_) * phi_pf * scalar_product(stress_term_du, z_vgrad) 
          + 2.0 * du_pf * (1.0 - constant_k_) * scalar_product(stress_term, z_vgrad) * phi_pf 
          + 2.0 * pf * (1.0 - constant_k_) * scalar_product(stress_term, z_vgrad) * du_pf) * state_fe_values.JxW(q_point);

          local_vector(i) += 0.0 * scale * (2.0 * pf * (1.0 - constant_k_) * scalar_product(stress_term_LinU, du_vgrad) * z_pf 
          + 2.0 * du_pf * (1.0 - constant_k_) * scalar_product(stress_term_LinU, grad_u) * z_pf 
          + 2.0 * scalar_product(stress_term_du, stress_term) * z_pf * phi_pf) * state_fe_values.JxW(q_point);
        }
      }
    }
  }

  void ElementEquation_Q(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double , double)
  {
    assert(this->problem_type_ == "gradient");
  }

  void ElementEquation_QT(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
    assert(this->problem_type_ == "tangent");
  }

  void ElementEquation_QTT(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
    assert(this->problem_type_ == "hessian");
  }

  void ElementEquation_QU(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
    assert(this->problem_type_ == "adjoint_hessian");
  }
  void ElementEquation_UQ(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
    assert(this->problem_type_ == "hessian");
  }
  void ElementEquation_QQ(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
    assert(this->problem_type_ == "hessian");
  }


  void BoundaryEquation(const FDC<DH, VECTOR, dealdim> &fdc, dealii::Vector<double> &local_vector, double scale, double)
  {
    assert(this->problem_type_ == "state");


    const auto &state_fe_face_values = fdc.GetFEFaceValuesState();
    unsigned int n_dofs_per_element = fdc.GetNDoFsPerElement();
    unsigned int n_q_points = fdc.GetNQPoints();
    unsigned int color = fdc.GetBoundaryIndicator();

    if (color == 3)
    {
      qvalues_.resize(n_q_points, Vector<double>(1));
      fdc.GetFaceValuesControl("control", qvalues_);

      const FEValuesExtractors::Vector displacement(0);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          const Tensor<1, 2> phi_i_v = state_fe_face_values[displacement].value(i, q_point);

          local_vector(i) -= scale * qvalues_[q_point][0] * phi_i_v[1] * state_fe_face_values.JxW(q_point);
        }
      }
    }
  }

  void BoundaryMatrix(const FDC<DH, VECTOR, dealdim> &, dealii::FullMatrix<double> &, double, double)
  {
  }

  void BoundaryRightHandSide(const FDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double)
  {
    assert(this->problem_type_ == "state");
  }

  void
  BoundaryEquation_Q(const FDC<DH, VECTOR, dealdim> &fdc, dealii::Vector<double> &local_vector, double scale, double)
  {
    assert(this->problem_type_ == "gradient");

    const auto &control_fe_face_values = fdc.GetFEFaceValuesControl();
    unsigned int n_dofs_per_element = fdc.GetNDoFsPerElement();
    unsigned int n_q_points = fdc.GetNQPoints();
    unsigned int color = fdc.GetBoundaryIndicator();

    zvalues_.resize(n_q_points, Vector<double>(3));
    fdc.GetFaceValuesState("adjoint", zvalues_);

    if (color == 3)
    {
      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          local_vector(i) += scale * (-zvalues_[q_point][1] * control_fe_face_values.shape_value(i, q_point)) 
          * control_fe_face_values.JxW(q_point);
        }
      }
    }
  }

  void BoundaryEquation_QT(const FDC<DH, VECTOR, dealdim> &fdc, dealii::Vector<double> &local_vector, double scale, double)
  {
    assert(this->problem_type_ == "tangent");

    const auto &state_fe_face_values = fdc.GetFEFaceValuesState();
    unsigned int n_dofs_per_element = fdc.GetNDoFsPerElement();
    unsigned int n_q_points = fdc.GetNQPoints();
    unsigned int color = fdc.GetBoundaryIndicator();

    if (color == 3)
    {
      dqvalues_.resize(n_q_points, Vector<double>(1));

      fdc.GetFaceValuesControl("dq", dqvalues_);

      const FEValuesExtractors::Vector displacement(0);

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          const Tensor<1, 2> phi_i_v = state_fe_face_values[displacement].value(i, q_point);

          local_vector(i) -= scale * dqvalues_[q_point][0] * phi_i_v[1] * state_fe_face_values.JxW(q_point);
        }
      }
    }
  }

  void BoundaryEquation_QTT(const FDC<DH, VECTOR, dealdim> &fdc, dealii::Vector<double> &local_vector, double scale, double)
  {
    assert(this->problem_type_ == "hessian");
    const auto &control_fe_face_values = fdc.GetFEFaceValuesControl();
    unsigned int n_dofs_per_element = fdc.GetNDoFsPerElement();
    unsigned int n_q_points = fdc.GetNQPoints();
    unsigned int color = fdc.GetBoundaryIndicator();

    dzvalues_.resize(n_q_points, Vector<double>(3));

    fdc.GetFaceValuesState("adjoint_hessian", dzvalues_);

    if (color == 3)
    {

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          local_vector(i) += scale * (-dzvalues_[q_point][1] * control_fe_face_values.shape_value(i, q_point)) * control_fe_face_values.JxW(q_point);
        }
      }
    }
  }

  void BoundaryEquation_U(const FDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
  }

  void BoundaryEquation_UT(const FDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
  }

  void BoundaryEquation_UTT(const FDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
  }

  void BoundaryEquation_UU(const FDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
  }

  void
  BoundaryEquation_QU(const FDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
  }

  void
  BoundaryEquation_UQ(const FDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
  }

  void
  BoundaryEquation_QQ(const FDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double, double)
  {
  }

  void ElementRightHandSide(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double)
  {
  }

  void ElementTimeEquation(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double)
  {
  }

  void ElementTimeMatrix(const EDC<DH, VECTOR, dealdim> &, FullMatrix<double> &)
  {
  }
  void ElementTimeEquation_U(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double)
  {
  }

  void ElementTimeEquation_UT(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double)
  {
  }

  void ElementTimeEquation_UTT(const EDC<DH, VECTOR, dealdim> &, dealii::Vector<double> &, double)
  {
  }

  void ElementTimeEquation_UU(const EDC<DH, VECTOR, dealdim> &edc, dealii::Vector<double> &local_vector, double scale)
  {
  }


  void ControlElementEquation(const EDC<DH, VECTOR, dealdim> &edc, dealii::Vector<double> &local_vector, double scale)
	{
    assert((this->problem_type_ == "gradient") || (this->problem_type_ == "hessian"));

    const DOpEWrapper::FEValues<dealdim> &control_fe_values = edc.GetFEValuesControl();
    unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
    unsigned int n_q_points = edc.GetNQPoints();

    funcgradvalues_.resize(n_q_points, Vector<double>(1));
    edc.GetValuesControl("last_newton_solution", funcgradvalues_);

    const FEValuesExtractors::Scalar controls(0);

    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
    {
      for (unsigned int i = 0; i < n_dofs_per_element; i++)
      {
        local_vector(i) += scale * (funcgradvalues_[q_point][0] * control_fe_values[controls].value(i, q_point)) * control_fe_values.JxW(q_point);
      }
    }
  }

  void ControlElementMatrix(const EDC<DH, VECTOR, dealdim> &edc,
                            FullMatrix<double> &local_matrix, double scale)
  {
    const DOpEWrapper::FEValues<dealdim> &control_fe_values = edc.GetFEValuesControl();
    unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
    unsigned int n_q_points = edc.GetNQPoints();

    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
    {
      for (unsigned int i = 0; i < n_dofs_per_element; i++)
      {
        for (unsigned int j = 0; j < n_dofs_per_element; j++)
        {
          local_matrix(i, j) += scale * control_fe_values.shape_value(i, q_point) 
          * control_fe_values.shape_value(j, q_point) * control_fe_values.JxW(q_point);
        }
      }
    }
  }

  void
  ControlBoundaryEquation(const FDC<DH, VECTOR, dealdim> &fdc, dealii::Vector<double> &local_vector, double scale)
  {
    assert((this->problem_type_ == "gradient") || (this->problem_type_ == "hessian"));

    const auto &control_fe_values = fdc.GetFEFaceValuesControl();
    unsigned int n_dofs_per_element = fdc.GetNDoFsPerElement();
    unsigned int n_q_points = fdc.GetNQPoints();
    unsigned int color = fdc.GetBoundaryIndicator();

    funcgradvalues_.resize(n_q_points, Vector<double>(1));
    fdc.GetFaceValuesControl("last_newton_solution", funcgradvalues_);

    const FEValuesExtractors::Scalar controls(0);

    if (color == 3)
    {
      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          local_vector(i) += scale * (funcgradvalues_[q_point][0] 
          * control_fe_values[controls].value(i, q_point)) * control_fe_values.JxW(q_point);
        }
      }
    }
  }

  void
  ControlBoundaryMatrix(const FDC<DH, VECTOR, dealdim> &fdc, FullMatrix<double> &local_matrix, double scale)
  {
    const auto &control_fe_values = fdc.GetFEFaceValuesControl();
    unsigned int n_dofs_per_element = fdc.GetNDoFsPerElement();
    unsigned int n_q_points = fdc.GetNQPoints();
    unsigned int color = fdc.GetBoundaryIndicator();

    if (color == 3)
    {
      for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
      {
        for (unsigned int i = 0; i < n_dofs_per_element; i++)
        {
          for (unsigned int j = 0; j < n_dofs_per_element; j++)
          {
            local_matrix(i, j) += scale *
                                  control_fe_values.shape_value(i, q_point) * control_fe_values.shape_value(j, q_point) * control_fe_values.JxW(q_point);
          }
        }
      }
    }
  }


  UpdateFlags GetUpdateFlags() const
  {
    return update_values | update_gradients | update_quadrature_points;
  }

  UpdateFlags GetFaceUpdateFlags() const
  {
    return update_values | update_gradients | update_normal_vectors | update_quadrature_points;
  }

  unsigned int
  GetControlNBlocks() const
  {
    return 1;
  }

  unsigned int
  GetStateNBlocks() const
  {
    return 2;
  }

  std::vector<unsigned int> &
  GetControlBlockComponent()
  {
    return control_block_component_;
  }
  const std::vector<unsigned int> &
  GetControlBlockComponent() const
  {
    return control_block_component_;
  }
  std::vector<unsigned int> &
  GetStateBlockComponent()
  {
    return state_block_component_;
  }
  const std::vector<unsigned int> &
  GetStateBlockComponent() const
  {
    return state_block_component_;
  }

private:
  std::vector<Vector<double>> funcgradvalues_;

  std::vector<Vector<double>> qvalues_;
  std::vector<Vector<double>> dqvalues_;

  std::vector<Vector<double>> uvalues_;
  std::vector<std::vector<Tensor<1, dealdim>>> ugrads_;
  std::vector<Vector<double>> last_timestep_uvalues_;
  std::vector<Vector<double>> next_timestep_uvalues_;

  std::vector<Vector<double>> zvalues_;
  std::vector<std::vector<Tensor<1, dealdim>>> zgrads_;
  std::vector<Vector<double>> last_timestep_zvalues_;

  std::vector<Vector<double>> duvalues_;
  std::vector<std::vector<Tensor<1, dealdim>>> dugrads_;
  std::vector<Vector<double>> last_timestep_duvalues_;

  std::vector<Vector<double>> dzvalues_;
  std::vector<std::vector<Tensor<1, dealdim>>> dzgrads_;
  std::vector<Vector<double>> last_timestep_dzvalues_;

  std::vector<std::vector<Tensor<1, dealdim>>> ufacegrads_;

  std::vector<unsigned int> state_block_component_;
  std::vector<unsigned int> control_block_component_;

  double constant_k_,
      alpha_eps_,
      G_c_,
      E_modulus_,
      poisson_ratio_nu_,
      lame_coefficient_mu_,
      lame_coefficient_lambda_,
      eta_, n_time_intervals_, end_time_,
      gamma_penal_, alpha_tik_,
      force_increment_;

  inline void Init_Displacement_Gradient(Tensor<2, 2> &T,
                                         const std::vector<Tensor<1, 2>> &u)
  {
    T.clear();
    T[0][0] = u[0][0];
    T[0][1] = 0.5 * u[0][1] + 0.5 * u[1][0];
    T[1][0] = T[0][1];
    T[1][1] = u[1][1];
  }
  
  inline void Init_Displacement(Tensor<1, 2> &T,
                                const Vector<double> &u)
  {
    T.clear();
    T[0] = u(0);
    T[1] = u(1);
  }

  inline void Init_Phasefield_Gradient(Tensor<1, 2> &T,
                                       const std::vector<Tensor<1, 2>> &u)
  {
    T.clear();
    T[0] = u[2][0];
    T[1] = u[2][1];
  }

  inline void Symmetrize(Tensor<2, 2> &T)
  {
    T[0][1] += T[1][0];
    T[0][1] *= 0.5;
    T[1][0] = T[0][1];
  }

  inline void Apply_Elasticity_Tensor(Tensor<2, 2> &T)
  {

    Tensor<2, 2> Identity;
    Identity[0][0] = 1.0;
    Identity[0][1] = 0.0;
    Identity[1][0] = 0.0;
    Identity[1][1] = 1.0;

    T = lame_coefficient_lambda_ * trace(T) * Identity +
        2.0 * lame_coefficient_mu_ * T;
  }
};
#endif
