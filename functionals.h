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

#ifndef LOCALFunctionalS_
#define LOCALFunctionalS_

#include <interfaces/pdeinterface.h>

using namespace dealii;
using namespace DOpE;


#if DEAL_II_VERSION_GTE(9,3,0)
template<
  template<bool DH, typename VECTOR, int dealdim> class EDC,
  template<bool DH, typename VECTOR, int dealdim> class FDC,
  bool DH, typename VECTOR, int dopedim, int dealdim>
class MinimalPhaseFieldValue : public FunctionalInterface<EDC, FDC, DH, VECTOR, dealdim>
#else
template<
  template<template<int, int> class DH, typename VECTOR, int dealdim> class EDC,
  template<template<int, int> class DH, typename VECTOR, int dealdim> class FDC,
  template<int, int> class DH, typename VECTOR, int dopedim, int dealdim>
class MinimalPhaseFieldValue : public FunctionalInterface<EDC, FDC, DH, VECTOR, dealdim>
#endif
{
    public:

        bool NeedTime() const
        {
            return true;
        }

        double ElementValue(const EDC<DH, VECTOR, dealdim>& edc)
        {
            unsigned int n_q_points = edc.GetNQPoints();
            double min_value_pf = 1.0e+16;

            std::vector<Vector<double> > uvalues;
            uvalues.resize(n_q_points, Vector<double>(3));
            edc.GetValuesState("state", uvalues);

            for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
            {
            double value_pf = uvalues[q_point](2);
            min_value_pf = value_pf;
            }
            return min_value_pf;
        }

        UpdateFlags GetUpdateFlags() const
        {
            return update_values | update_quadrature_points;
        }

        std::string GetType() const
        {
            return "domain timelocal";
        }
        std::string GetName() const
        {
            return "Phase-field summed up";
        }
};


#if DEAL_II_VERSION_GTE(9,3,0)
template<
  template<bool DH, typename VECTOR, int dealdim> class EDC,
  template<bool DH, typename VECTOR, int dealdim> class FDC,
  bool DH, typename VECTOR, int dopedim, int dealdim>
class Term_phi_phi_d : public FunctionalInterface<EDC, FDC, DH,VECTOR, dealdim>
#else
template<
  template<template<int, int> class DH, typename VECTOR, int dealdim> class EDC,
  template<template<int, int> class DH, typename VECTOR, int dealdim> class FDC,
  template<int, int> class DH, typename VECTOR, int dopedim, int dealdim>
class Term_phi_phi_d : public FunctionalInterface<EDC, FDC, DH,VECTOR, dealdim>
#endif
{
    public:

        bool NeedTime() const
        {
            if(this->GetTime() > 0.)
            return true;
            return false;
        }

        void update_Crack(double x_min, double x_max, double y_min, double y_max)
        {
            crack_low_x_ = x_min;
            crack_high_x_ = x_max;
            crack_low_y_ = y_min;
            crack_high_y_ = y_max;
        }

        double ElementValue(const EDC<DH, VECTOR, dealdim>& edc)
        {
            const DOpEWrapper::FEValues<dealdim> &state_fe_values = edc.GetFEValuesState();
            unsigned int n_q_points = edc.GetNQPoints();

            std::vector<Vector<double> > uvalues_;
            std::vector<Vector<double> > phi_d_values;
            uvalues_.resize(n_q_points, Vector<double>(3));
            edc.GetValuesState("state", uvalues_);
            phi_d_values.resize(n_q_points, Vector<double>(3));

            double pf = 0.;
            double ret = 0.;

            for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
            {
            
                pf = uvalues_[q_point](2);

                phi_d_values[q_point][2] =1.0;
                if(state_fe_values.quadrature_point(q_point)[1]>=crack_low_y_ 
                &&state_fe_values.quadrature_point(q_point)[1]<=crack_high_y_
                &&state_fe_values.quadrature_point(q_point)[0]<=crack_high_x_
                &&state_fe_values.quadrature_point(q_point)[0]>= crack_low_x_)
                {
                    phi_d_values[q_point][2]=0;
                }

                ret += 0.5 * (pf - phi_d_values[q_point][2]) 
                           * (pf-phi_d_values[q_point][2])
                           * state_fe_values.JxW(q_point);
            }
            return ret;
        }

    UpdateFlags GetUpdateFlags() const
    {
        return update_values | update_quadrature_points;
    }

    std::string GetType() const
    {
        return "domain timedistributed";
    }

    std::string GetName() const
    {
        return "Cost_Term_phi_d";
    }

    private:
        double gamma_penal_, crack_low_x_, crack_high_x_, crack_low_y_, crack_high_y_;
};
#endif
