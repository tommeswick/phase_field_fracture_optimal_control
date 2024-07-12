/**
*
* Copyright (C) 2012-2014 by the DOpElib authors
*
* This file is part of DOpElib
*
* DOpElib is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* DOpElib is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT included in this distribution
* for further information on this license.
*
**/

#include <wrapper/function_wrapper.h>

using namespace dealii;



/******************************************************/


class InitialData : public DOpEWrapper::Function<2>
{
  public:
    InitialData(const double cell_diameter) :
        DOpEWrapper::Function<2>(3)
    {
      cell_diameter_ = cell_diameter;
    }
    virtual double
    value(const Point<2> &p, const unsigned int component = 0) const;
    virtual void
    vector_value(const Point<2> &p, Vector<double> &value) const;

  private:
    double cell_diameter_;

};

/******************************************************/

double
InitialData::value(const Point<2> &p, const unsigned int component) const
{
  // TODO: should be the local element size
  double min_dia_ = cell_diameter_;
  // Sanity check of the diameter
  //std::cout << "Cell diameter: " << min_dia_ << std::endl;

   if (component == 2)
      {
	if (((p(0) >= -0.1  - min_dia_) && (p(0) <= 0.1 + min_dia_ )) &&
	    ((p(1) >=  0.0 - min_dia_) && (p(1) <=  0.0 + min_dia_))
	    )
	  return 1.0;
	else
	  return 1.0;
      }

    return 0.0;

}

/******************************************************/

void
InitialData::vector_value(const Point<2> &p, Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = InitialData::value(p, c);
}

/******************************************************/

class BoundaryParabel : public DOpEWrapper::Function<2>
{
public:
  BoundaryParabel (ParameterReader &param_reader,
		   const double cell_diameter) : DOpEWrapper::Function<2>(3)
  {
    param_reader.SetSubsection("My functions parameters");
    displacement_step_size_a_ = param_reader.get_double ("displacement_step_size_a");
    displacement_step_size_b_ = param_reader.get_double ("displacement_step_size_b");

    cell_diameter_ = cell_diameter;
  }

  virtual double value (const Point<2>   &p,
			const unsigned int  component = 0) const;

  virtual void vector_value (const Point<2> &p,
			     Vector<double>   &value) const;

  static void declare_params(ParameterReader &param_reader)
  {
    param_reader.SetSubsection("My functions parameters");
    param_reader.declare_entry("displacement_step_size_a", "0.0",
			       Patterns::Double(0));
    param_reader.declare_entry("displacement_step_size_b", "0.0",
			       Patterns::Double(0));
  }

 void SetTime(double t) const
 {
   localtime=t;
 }

private:
 double cell_diameter_, displacement_step_size_a_, displacement_step_size_b_;
  mutable double localtime;



};

/******************************************************/

double
BoundaryParabel::value (const Point<2>  &p,
			    const unsigned int component) const
{
  Assert (component < this->n_components,
	  ExcIndexRange (component, 0, this->n_components));

  //std::cout << localtime << std::endl;
  if (component == 0)
    {
      if (localtime < 200)
	{
	  return ( ((p(1) == 1.0) && (p(0) <= 1.0) && (p(0) >= -1.0))
		   ?
		   (1.0) * localtime * displacement_step_size_a_ : 0 );
	}
      else
	{
	  return ( ((p(1) == 1.0) && (p(0) <= 1.0) && (p(0) >= -1.0))
		   ?
		   (1.0) * localtime * displacement_step_size_b_ : 0 );
	}

    }

  if (component == 2)
    {
      return 1.0;
    }

  return 0;
}

/******************************************************/

void
BoundaryParabel::vector_value (const Point<2> &p,
				   Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values (c) = BoundaryParabel::value (p, c);
}
