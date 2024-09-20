/** Copyright (C) 2024 by Denis Khimin, Marc C. Steinbach, Thomas Wick   
 *
 *  This file is part of the phase_field_fracture_optimal_control library.
 *
 *  The phase_field_fracture_optimal_control library is free software; 
 *  you can use it, redistribute it, and/or modify it
 *  under the terms of the GNU Lesser General
 *  Public License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *  The full text of the license can be found in the file LICENSE. 
**/

#ifndef LOCALFunctional_
#define LOCALFunctional_

#include <interfaces/functionalinterface.h>

using namespace dealii;
using namespace DOpE;

#if DEAL_II_VERSION_GTE(9,3,0)
template<
  template<bool DH, typename VECTOR, int dealdim> class EDC,
  template<bool DH, typename VECTOR, int dealdim> class FDC,
  bool DH, typename VECTOR, int dopedim, int dealdim>
class LocalFunctional : public FunctionalInterface<EDC, FDC, DH, VECTOR,
  dopedim, dealdim>
#else
template<
  template<template<int, int> class DH, typename VECTOR, int dealdim> class EDC,
  template<template<int, int> class DH, typename VECTOR, int dealdim> class FDC,
  template<int, int> class DH, typename VECTOR, int dopedim, int dealdim>
class LocalFunctional : public FunctionalInterface<EDC, FDC, DH, VECTOR,
  dopedim, dealdim>
#endif
{
public:

  static void declare_params(ParameterReader &param_reader)
  {
    param_reader.SetSubsection("Local Functional parameters");
    param_reader.declare_entry("alpha_tikhonov", "0.0", Patterns::Double(0));
    param_reader.declare_entry("control_constant", "0.0", Patterns::Double(0));
    param_reader.declare_entry("ud_values_constant", "0.0", Patterns::Double(0));
    param_reader.declare_entry("ud_values_scale", "0.0", Patterns::Double(0));

    param_reader.SetSubsection("Force parameters");
    param_reader.declare_entry("force_increment", "0.0", Patterns::Double(0));

    param_reader.SetSubsection("Crack parameters");
    param_reader.declare_entry("crack_low_x", "0.0", Patterns::Double(0));
    param_reader.declare_entry("crack_high_x", "0.0", Patterns::Double(0));
    param_reader.declare_entry("crack_low_y", "0.0", Patterns::Double(0));
    param_reader.declare_entry("crack_high_y", "0.0", Patterns::Double(0));
  }

  void update_Crack(double x_min, double x_max, double y_min, double y_max)
  {
    crack_low_x_ = x_min;
    crack_high_x_ = x_max;
    crack_low_y_ = y_min;
    crack_high_y_ = y_max;
  }

  LocalFunctional(ParameterReader &param_reader)
  {
		param_reader.SetSubsection("Local Functional parameters");
		alpha_tikhonov_ = param_reader.get_double("alpha_tikhonov");
		control_constant_ = param_reader.get_double("control_constant");
		ud_values_constant_ = param_reader.get_double("ud_values_constant");
		ud_values_scale_ = param_reader.get_double("ud_values_scale");

		param_reader.SetSubsection("Crack parameters");
		crack_low_x_ = param_reader.get_double("crack_low_x");
		crack_high_x_ = param_reader.get_double("crack_high_x");
		crack_low_y_ = param_reader.get_double("crack_low_y");
		crack_high_y_ = param_reader.get_double("crack_high_y");
  }

  bool NeedTime() const
  {
    if (this->GetTime() > 0.)
      return true;
    return false;
  }

  double ElementValue(const EDC<DH, VECTOR, dealdim>& edc)
  {
    const DOpEWrapper::FEValues<dealdim> & state_fe_values = edc.GetFEValuesState();
    unsigned int n_q_points = edc.GetNQPoints();

    double ret = 0.;

    phi_d_values.resize(n_q_points, Vector<double>(3));
    uvalues_.resize(n_q_points, Vector<double>(3));
    edc.GetValuesState("state", uvalues_);

    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
		{
				double pf = uvalues_[q_point](2);

				phi_d_values[q_point][2] =1.0;

				if(state_fe_values.quadrature_point(q_point)[1]>=crack_low_y_
				 &&state_fe_values.quadrature_point(q_point)[1]<=crack_high_y_
				 &&state_fe_values.quadrature_point(q_point)[0]<=crack_high_x_ 
			   &&state_fe_values.quadrature_point(q_point)[0]>= crack_low_x_)
				{
	      	phi_d_values[q_point][2]=0;
	    	}
			ret += 0.5 * (pf-phi_d_values[q_point][2])* (pf-phi_d_values[q_point][2])* state_fe_values.JxW(q_point);
		}
  	return ret;
	}
  void ElementValue_U(const EDC<DH, VECTOR, dealdim>& edc, dealii::Vector<double> &local_vector, double scale)
  {
    const DOpEWrapper::FEValues<dealdim> & state_fe_values = edc.GetFEValuesState();
    unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
    unsigned int n_q_points = edc.GetNQPoints();

    uvalues_.resize(n_q_points, Vector<double>(3));
    ugrads_.resize(n_q_points, std::vector<Tensor<1, 2> >(3));
    phi_d_values.resize(n_q_points, Vector<double>(3));
    edc.GetValuesState("state", uvalues_);
    edc.GetGradsState("state", ugrads_);

    const FEValuesExtractors::Vector displacement(0);
    const FEValuesExtractors::Scalar phase(2);

		for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
		{
	
			double pf = uvalues_[q_point](2);

			phi_d_values[q_point][2] = 1;
			if(state_fe_values.quadrature_point(q_point)[1]>=crack_low_y_
				 &&state_fe_values.quadrature_point(q_point)[1]<=crack_high_y_
				 &&state_fe_values.quadrature_point(q_point)[0]<=crack_high_x_ 
			   &&state_fe_values.quadrature_point(q_point)[0]>= crack_low_x_)
				{
	      	phi_d_values[q_point][2]=0;
	    	}

			for (unsigned int i = 0; i < n_dofs_per_element; i++)
	  	{
	    	double phi_i_pf = state_fe_values[phase].value(i, q_point);
	    	local_vector(i) += scale*0.5*((pf-phi_d_values[q_point][2])*phi_i_pf)*state_fe_values.JxW(q_point);
	  	}
		}
  }
 
  void ElementValue_UU(const EDC<DH, VECTOR, dealdim>& edc, dealii::Vector<double> &local_vector, double scale)
  {
    const DOpEWrapper::FEValues<dealdim> & state_fe_values = edc.GetFEValuesState();
    unsigned int n_dofs_per_element = edc.GetNDoFsPerElement();
    unsigned int n_q_points = edc.GetNQPoints();

    uvalues_.resize(n_q_points, Vector<double>(3));
    duvalues_.resize(n_q_points, Vector<double>(3));
    edc.GetValuesState("state", uvalues_);
    edc.GetValuesState("tangent", duvalues_);

    const FEValuesExtractors::Vector displacement(0);
    const FEValuesExtractors::Scalar phase(2);

    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
    {
			double du_pf = duvalues_[q_point](2);
			for (unsigned int i = 0; i < n_dofs_per_element; i++)
	  	{
	    	double phi_pf = state_fe_values[phase].value(i, q_point);
	    	local_vector(i) += scale*0.5 * du_pf *phi_pf *state_fe_values.JxW(q_point);
	  	}
    }
  }

  void ElementValue_Q(const EDC<DH, VECTOR, dealdim>&, dealii::Vector<double> &, double)
  {
  }

  void ElementValue_QU(const EDC<DH, VECTOR, dealdim>&, dealii::Vector<double> &, double)
  {
  }

  void
  ElementValue_UQ(const EDC<DH, VECTOR, dealdim>&, dealii::Vector<double> &, double)
  {
  }

  void
  ElementValue_QQ(const EDC<DH, VECTOR, dealdim>&, dealii::Vector<double> &, double)
  {
  }
  
  double BoundaryValue(const FDC<DH, VECTOR, dealdim>& fdc)
  {
    const auto & state_fe_values = fdc.GetFEFaceValuesControl();
    unsigned int n_q_points = fdc.GetNQPoints();

    double ret = 0.;

    qvalues_.resize(n_q_points, Vector<double>(1));
    fdc.GetFaceValuesControl("control", qvalues_);

    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
    {
			ret += 0.5 * alpha_tikhonov_ * (qvalues_[q_point][0]  - control_constant_)
			*(qvalues_[q_point][0] - control_constant_) * state_fe_values.JxW(q_point);
    }
    return ret;
  }

  void BoundaryValue_Q(const FDC<DH, VECTOR, dealdim>& fdc, dealii::Vector<double> &local_vector, double scale)
  {
    const auto & state_fe_values =fdc.GetFEFaceValuesControl();
    unsigned int n_dofs_per_element = fdc.GetNDoFsPerElement();
    unsigned int n_q_points = fdc.GetNQPoints();

    qvalues_.resize(n_q_points, Vector<double>(1));
    fdc.GetFaceValuesControl("control", qvalues_);

    const FEValuesExtractors::Scalar control_ext(0);
    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
    {
			for (unsigned int i = 0; i < n_dofs_per_element; i++)
	  	{
	    	local_vector(i) += scale * alpha_tikhonov_ * (qvalues_[q_point][0] - control_constant_)
	      * state_fe_values[control_ext].value(i, q_point) * state_fe_values.JxW(q_point);
	  	}
		}
  }

  void BoundaryValue_QQ(const FDC<DH, VECTOR, dealdim>& fdc, dealii::Vector<double> &local_vector, double scale)
  {
    const auto & state_fe_values = fdc.GetFEFaceValuesControl();
    unsigned int n_dofs_per_element = fdc.GetNDoFsPerElement();
    unsigned int n_q_points = fdc.GetNQPoints();

    dqvalues_.resize(n_q_points, Vector<double>(1));
    fdc.GetFaceValuesControl("dq", dqvalues_);

    const FEValuesExtractors::Scalar control_ext(0);
    for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
    {
			for (unsigned int i = 0; i < n_dofs_per_element; i++)
	  	{
	    	local_vector(i) += scale * alpha_tikhonov_ * dqvalues_[q_point][0] * state_fe_values[control_ext].value(i, q_point)
	      * state_fe_values.JxW(q_point);
	  	}
    }
  }

  void BoundaryValue_U(const FDC<DH, VECTOR, dealdim>&, dealii::Vector<double> &, double)
  {
  }

  void BoundaryValue_UU(const FDC<DH, VECTOR, dealdim>&, dealii::Vector<double> &, double)
  {
  }
  
  void BoundaryValue_QU(const FDC<DH, VECTOR, dealdim>&, dealii::Vector<double> &, double)
  {

  }
  void BoundaryValue_UQ(const FDC<DH, VECTOR, dealdim>&, dealii::Vector<double> &, double)
  {
  }

  UpdateFlags GetUpdateFlags() const
  {
    return update_quadrature_points | update_values |update_gradients;
  }

  UpdateFlags GetFaceUpdateFlags() const
  {
    return update_values | update_normal_vectors | update_quadrature_points | update_gradients;
  }

  std::string GetType() const
  {
    return "domain boundary timedistributed";
  }

  std::string GetName() const
  {
    return "Cost functional";
  }

private:
  std::vector<Vector<double> > uvalues_;
  std::vector<std::vector<Tensor<1, dealdim> > > ugrads_;
  std::vector<Vector<double> > phi_d_values;
  std::vector<Vector<double> > duvalues_;
  std::vector<std::vector<Tensor<1, dealdim> > > dugrads_;

  std::vector<Vector<double> > qvalues_;
  std::vector<Vector<double> > dqvalues_;

  double alpha_tikhonov_,
  control_constant_,
  ud_values_constant_,
  ud_values_scale_,
  time_i,
  crack_low_x_,
  crack_high_x_,
  crack_low_y_,
  crack_high_y_;

  inline void Init_Displacement(Tensor<1,2> & T, const Vector<double> & u)
  {
    T.clear();
    T[0] = u(0);
    T[1] = u(1);
  }
  inline void Init_Phasefield_Gradient(Tensor<1,2> & T, const std::vector<Tensor<1,2> >& u)
  {
    T.clear();
    T[0] = u[2][0];
    T[1] = u[2][1];
  }
};
#endif
