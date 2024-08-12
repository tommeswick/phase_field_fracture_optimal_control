/** Authors: Denis Khimin, Marc C. Steinbach, Thomas Wick 
 *  Leibniz University Hannover  
 *  Institute of Applied Mathematics
 *  Date: Jul 12, 2024 
 *                                                                
 *  Copyright (C) 2024 by Denis Khimin, Marc C. Steinbach, Thomas Wick                    
 *                                                             
 *  Journal Reference: 
 *  authors: Denis Khimin, Marc. C. Steinbach, Thomas Wick
 *  title:   Space-time formulation, discretization, and computational performance studies 
 *           for phase-field fracture optimal control problems (Sec. 5.1.1)
 *  journal: Journal of Computational Physics (JCP), Vol. 470, 2022, 111554
 *  link:    https://doi.org/10.1016/j.jcp.2022.111554
**/ 

/** Library information
 *  -------------------
 *  1. Install deal.II 9.5.1 via www.dealii.org
 *     Download: https://www.dealii.org/download.html
 *     Installation instructions: https://www.dealii.org/current/readme.html
 *
 *  2. Install DOpElib via http://www.dopelib.net
 *     Download: https://github.com/winnifried/dopelib
 *     Installation instructions: https://winnifried.github.io/dopelib/documentation.html
 *       in Chapter 2 of the *.pdf manual.
 *
 *  3. Please put the current code into some new folder on your machine.
 *     Follow instructions from DOpElib *.pdf manual in Chapter 4 
 *     (i.e, Section 4.4 Creating new examples)
 *     to set up all environment variables (for finding deal.II and DOpElib) correctly.
 *     Then: build, compile, run as described in Section 4.4 of the dopelib manual.
 *
 *  4. The results of this code (see local folder Results/ ) should then reproduce 
 *     Example 1 (Section 5.1.1) of Khimin et al., JCP, 2022. 
**/


// Include standard libraries
#include <iostream>
#include <fstream>

// Include deal.II libraries
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#if DEAL_II_VERSION_GTE(9,1,1)
#else
#include <deal.II/grid/tria_boundary_lib.h>
#endif
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

// Include DOpE libraries
#include <include/parameterreader.h>
#include <templates/directlinearsolver.h>
#include <templates/integrator.h>
#include <basic/mol_spacetimehandler.h>
#include <problemdata/simpledirichletdata.h>
#include <container/integratordatacontainer.h>
#include <templates/newtonsolver.h>
#include <interfaces/functionalinterface.h>
#include <problemdata/noconstraints.h>
#include <opt_algorithms/reducednewtonalgorithm.h>

// Include DOpE libraries for instationary problems
#include <reducedproblems/instatreducedproblem.h>
//#include <templates/instat_step_newtonsolver.h>
#include "instat_step_modified_newtonsolver.h"
#include <container/instatoptproblemcontainer.h>
#include <tsschemes/backward_euler_problem.h>
#include <tsschemes/shifted_crank_nicolson_problem.h>

// Include problem-specific headers from this local folder.
// These files contain the models and parameters that
// can be changed by the user depending on the given problem
// statement.
#include "localpde_eta.h"
#include "localfunctional.h"
#include "functionals.h"
#include "my_functions.h"

using namespace dealii;
using namespace DOpE;

// Dimension of the FEM forward problem
const static int DIM = 2;
// Dimension of the FEM control problem
const static int CDIM = 2;

// Macros
#if DEAL_II_VERSION_GTE(9,3,0)
#define DOFHANDLER false
#else
#define DOFHANDLER DoFHandler
#endif
#define FE FESystem
#define CDC ElementDataContainer
#define FDC FaceDataContainer

typedef QGauss<DIM> QUADRATURE;
typedef QGauss<DIM - 1> FACEQUADRATURE;
typedef BlockSparseMatrix<double> MATRIX;
typedef BlockSparsityPattern SPARSITYPATTERN;
typedef BlockVector<double> VECTOR;

// Define the functional interface
typedef FunctionalInterface<CDC, FDC, DOFHANDLER, VECTOR, CDIM, DIM> FUNC;

// Define the optimization problem container
typedef OptProblemContainer<
    LocalFunctional<CDC, FDC, DOFHANDLER, VECTOR, CDIM, DIM>, FUNC,
    LocalPDE<CDC, FDC, DOFHANDLER, VECTOR, DIM>,
    SimpleDirichletData<VECTOR, DIM>,
    NoConstraints<CDC, FDC, DOFHANDLER, VECTOR, CDIM, DIM>, SPARSITYPATTERN,
    VECTOR, CDIM, DIM> OP_BASE;

// Define the state problem, i.e. the PDE constraint of the optimization problem
// In our case this will be phase-field fracture
typedef StateProblem<OP_BASE, LocalPDE<CDC, FDC, DOFHANDLER, VECTOR, DIM>,
    SimpleDirichletData<VECTOR, DIM>, SPARSITYPATTERN, VECTOR, DIM> PROB;

// Define macros for time-stepping schemes
#define TSP BackwardEulerProblem
#define DTSP BackwardEulerProblem
typedef InstatOptProblemContainer<TSP, DTSP, FUNC,
    LocalFunctional<CDC, FDC, DOFHANDLER, VECTOR, CDIM, DIM>,
    LocalPDE<CDC, FDC, DOFHANDLER, VECTOR, DIM>,
    SimpleDirichletData<VECTOR, DIM>,
    NoConstraints<CDC, FDC, DOFHANDLER, VECTOR, DIM, DIM>, SPARSITYPATTERN,
    VECTOR, DIM, DIM> OP;
#undef TSP
#undef DTSP

// Define integrator and the NLP (nonlinear programming) solver
typedef IntegratorDataContainer<DOFHANDLER, QUADRATURE,
    FACEQUADRATURE, VECTOR, DIM> IDC;
typedef Integrator<IDC, VECTOR, double, DIM> INTEGRATOR;
typedef DirectLinearSolverWithMatrix<SPARSITYPATTERN, MATRIX, VECTOR> LINEARSOLVER;
typedef NewtonSolver<INTEGRATOR, LINEARSOLVER, VECTOR> CNLS;

// Define the Newton solver which is used for all PDE problems,
// i.e., this Newton solver treats the nonlinear forward PDE problem.
typedef InstatStepModifiedNewtonSolver<INTEGRATOR, LINEARSOLVER, VECTOR> NLS;

/** 
 * 
 * Define a type for the reduced Newton algorithm. Within the reduced approach, 
 * the original cost functional J(q,U), depending on the control q and the state U
 * is replaced with the reduced cost functional j(q) = J(q,S(q)), where S
 * is the solution operator for the PDE problem. Hence the reduced 
 * Newton solves j'(q) = 0.
 * For further details see Khimin et al., JCP, 2022, Section 3.2 and Algorithm 1.
 * 
**/
typedef ReducedNewtonAlgorithm<OP, VECTOR> RNA;
typedef InstatReducedProblem<CNLS, NLS, INTEGRATOR, INTEGRATOR, OP, VECTOR, DIM,
    DIM> RP;

/////////////////////////////////////////////////////////////////////////////////
// Declare mathematical, physical and numerical parameters, 
// which are read from the parameter file dope.prm
void
declare_params(ParameterReader &param_reader)
{
  param_reader.SetSubsection("main parameters");
  param_reader.declare_entry("initial_control", "0.0", Patterns::Double(0));
  param_reader.declare_entry("time_intervals", "0", Patterns::Integer(0));
  param_reader.declare_entry("end_time", "0", Patterns::Double(0));
  param_reader.declare_entry("refinement_level", "6", Patterns::Integer(0));

	param_reader.SetSubsection("Local Functional parameters");
	param_reader.declare_entry("alpha_tikhonov", "0.0", Patterns::Double(0));
	param_reader.declare_entry("control_constant", "0.0", Patterns::Double(0));

  param_reader.SetSubsection("Local PDE parameters");
  param_reader.declare_entry("constant_k", "0.0", Patterns::Double(0));
	param_reader.declare_entry("alpha_eps", "0.0", Patterns::Double(0));
	param_reader.declare_entry("gamma_penal", "0.0", Patterns::Double(0));
	param_reader.declare_entry("G_c", "0.0", Patterns::Double(0));
	param_reader.declare_entry("E_modulus", "0.0", Patterns::Double(0));
	param_reader.declare_entry("poisson_ratio_nu", "0.0", Patterns::Double(0));
	param_reader.declare_entry("eta", "0.0", Patterns::Double(0));
  /**
   * 
   * The meaning of the parameters is explained 
   * in the introduction and in chap 5 of Khimin et al., JCP, 2022.
   * 
  **/
}


/////////////////////////////////////////////////////////////////////////////////
// Main function
int main(int argc, char **argv)
{
  // Initialize MPI
  dealii::Utilities::MPI::MPI_InitFinalize mpi(argc, argv);

  // Read parameter file
  std::string paramfile = "dope.prm";
  if (argc == 2)
  {
    paramfile = argv[1];
  }
  else if (argc > 2)
  {
    std::cout << "Usage: " << argv[0] << " [ paramfile ] " << std::endl;
    return -1;
  }

  // Initialize parameter reader and declare parameters
  ParameterReader pr;
  RP::declare_params(pr);
  RNA::declare_params(pr);
  LocalPDE<CDC, FDC, DOFHANDLER, VECTOR, DIM>::declare_params(pr);
  LocalFunctional<CDC, FDC, DOFHANDLER, VECTOR, CDIM, DIM>::declare_params(pr);
  BoundaryParabel::declare_params(pr);
  declare_params(pr);
  pr.read_parameters(paramfile);

  // Read main parameters
  pr.SetSubsection("main parameters");
  const double initial_control = pr.get_double("initial_control");
  const unsigned int time_intervals = pr.get_integer("time_intervals");
  const double end_time = pr.get_double("end_time");
  const unsigned int refinement_level = pr.get_integer("refinement_level");

  // Setup triangulation (e.g., the mesh) and read input file
  Triangulation<DIM> triangulation;
  GridIn<DIM> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream input_file("unit_slit.inp");
  grid_in.read_ucd(input_file);

  // Refine grid uniformly
  triangulation.refine_global(refinement_level);
  const double cell_diameter = 
	GridTools::diameter(triangulation)/std::pow(2.0,refinement_level);
  
/** 
   * 
   * Define finite elements, i.e., we use Q1 finite elements in 2D
   * for the control and for the solution/state variable, 
   * which has two entries U = (u,phi) = (u_1,u_2,phi), 
   * see Sec. 2.1 (Khimin et al., JCP, 2022)
   * 
  **/
  FE<DIM> control_fe(FE_Q<CDIM>(1), 1);
  FE<DIM> state_fe(FE_Q<DIM>(1), 2, FE_Q<DIM>(1), 1);

  // Define quadrature formulas to evaluate FEM integrals
  QUADRATURE quadrature_formula(3);
  FACEQUADRATURE face_quadrature_formula(3);

  // Initialize integrator data container
  IDC idc(quadrature_formula, face_quadrature_formula);

  // Initialize local PDE
  LocalPDE<CDC, FDC, DOFHANDLER, VECTOR, DIM> LPDE(pr);

  // Read and set local functional parameters of 
  // the optimal control part
	pr.SetSubsection("Local Functional parameters");
	double alpha_tikhonov = pr.get_double("alpha_tikhonov");
	double control_constant = pr.get_double("control_constant");

  // Read and set local PDE parameters of the phase-field fracture forward
  // problem.
  pr.SetSubsection("Local PDE parameters");
  double gamma_penal = pr.get_double("gamma_penal");
	double constant_k = pr.get_double("constant_k");
	double alpha_eps = pr.get_double("alpha_eps");
	double G_c = pr.get_double("G_c");
	double E_modulus = pr.get_double("E_modulus");
	double poisson_ratio_nu = pr.get_double("poisson_ratio_nu");
	double eta = pr.get_double("eta");

  /**
   * 
   * Read and set the desired crack, which is given as a rectangle
   * consisting of all points (x,y) that satisfy
   * crack_low_x <= x <= crack_high_x
   * crack_low_y <= y <= crack_high_y
   * 
  **/
  pr.SetSubsection("Crack parameters");
  double crack_low_x = pr.get_double("crack_low_x");
	double crack_high_x = pr.get_double("crack_high_x");
	double crack_low_y = pr.get_double("crack_low_y");
	double crack_high_y = pr.get_double("crack_high_y");

  // Update gamma penalization parameter for local PDE
  LPDE.update_gamma(gamma_penal);

/**
   * 
   * Initialize local functional and minimal phase field value.
   * Localfunctional --> Cost functional J (sec. 3.1 of Khimin et al. JCP, 2022)
   * MinimalPhaseFieldValues --> Measures the propagation of the crack
   * Term_phi_phi_d --> Measures the tracking part of J
   * 
  **/
  LocalFunctional<CDC, FDC, DOFHANDLER, VECTOR, CDIM, DIM> LFunc(pr);
  MinimalPhaseFieldValue <CDC, FDC, DOFHANDLER, VECTOR, CDIM, DIM> MPFV;
  Term_phi_phi_d <CDC, FDC, DOFHANDLER, VECTOR, CDIM, DIM> FUNC4;

  // Define time intervals for the problem. We work in a space-time 
  // formulation and consequently, the time interval is inialized 
  // similar to the spatial triangulation.
  Triangulation<1> times;
  GridGenerator::subdivided_hyper_cube(times, time_intervals, 0, end_time);

  // Initialize Method of Lines handler
  // At this point, we have discretized the time domain with dG(0) and obtain
  // a backward Euler scheme. This results in a classical methods of line 
  // time discretization.
  MethodOfLines_SpaceTimeHandler<FE, DOFHANDLER, SPARSITYPATTERN, VECTOR, DIM,
      DIM> DOFH(triangulation, control_fe, state_fe, times,
      DOpEtypes::VectorAction::stationary);

  // Initialize constraints (for examples and more information,
  // please refer to the DOpElib documentation)
  NoConstraints<CDC, FDC, DOFHANDLER, VECTOR, DIM, DIM> Constraints;

  // Initialize optimization problem
  OP P(LFunc, LPDE, Constraints, DOFH);
  P.AddFunctional(&MPFV);
  P.AddFunctional(&FUNC4);

  /**
   * 
   * Set boundary conditions and specify the boundary where the
   * control should act. Here it acts on the top boundary with id = 3.
   * The comp_mask specifies to which component the boundary conditions
   * should be applied: u_1,u_2, or phi.
   * 
  **/ 
  std::vector<bool> comp_mask(3);
  comp_mask[0] = true;
  comp_mask[1] = true;
  comp_mask[2] = false;

  // Initialize zero function for Dirichlet boundary data as DD1
  DOpEWrapper::ZeroFunction<DIM> zf(3);
  SimpleDirichletData<VECTOR, DIM> DD1(zf);

  // Initialize boundary parabel and Dirichlet boundary data as DD2
  BoundaryParabel boundary_parabel(pr,cell_diameter);
  SimpleDirichletData<VECTOR, DIM> DD2(boundary_parabel);

  // Apply homogeneous boundary values zu u_1 and u_2 on the
  // bottom boundary with id = 2
  comp_mask[0] = true;
  comp_mask[1] = true;
  P.SetDirichletBoundaryColors(2, comp_mask, &DD1);

  // The control and the corresponding equations act only on 
  // the top boundary with id = 3
  P.SetBoundaryEquationColors(3);
  P.SetBoundaryFunctionalColors(3);
  P.SetControlBoundaryEquationColors(3); 

  // Set initial values
  InitialData initial_data(cell_diameter);
  P.SetInitialValues(&initial_data);


  // Initialize solver and algorithm, where Alg is given in Sec. 4.5.
  RP solver(&P, DOpEtypes::VectorStorageType::fullmem, pr, idc);
  RNA Alg(&P, &solver, pr);
  Alg.ReInit();

  // Initialize control vector and set it to the initial control value
  ControlVector<VECTOR> q(&DOFH, DOpEtypes::VectorStorageType::fullmem,pr);
  q = initial_control;

  // Since FUNC4 evaluates the tracking term which involves the desired
  // phase-field, we have to update the parameters defining the desired crack
  // in this function
	FUNC4.update_Crack(crack_low_x, crack_high_x, crack_low_y, crack_high_y);

  // Output
	std::stringstream out;
	out << "\n" << "Initial control: " << initial_control << "\n";
	out << "Time intervals: " << time_intervals << "\n";
	out << "Cell diameter: " << cell_diameter << "\n";
	out << "Tikhonov parameter: " << alpha_tikhonov << "\n";
	out << "Control constant: " << control_constant << "\n";
	out << "Kappa: " << constant_k << "\n";
	out << "Epsilon: " << alpha_eps << "\n";
  out << "Gamma: " << gamma_penal << "\n";
	out << "Fracture Toughness: " << G_c << "\n";
	out << "Young's modulus: " << E_modulus << "\n";
  out << "Poisson ratio: " << poisson_ratio_nu << "\n";
	out << "Viscous parameter: " << eta << "\n";
	out<<"INITIAL DESIRED CRACK varphi_d := ["<<crack_low_x<<","<<
	crack_high_x<<"] x ["<<crack_low_y<<","<<crack_high_y<<"]."<<"\n";

/**
   *
   * In Sec. 5.5.1 (Khimin et al., JCP, 2022) our approch is to repeatedly solve the optimal control problem
   * defined in P. In each iteration we update the desired phase-field, 
   * i.e., we prolong it to the left, and 
   * solve the optimal control problem again with the new desired phase-field. 
   * In Sec. 5.1.1 the left boundary of the desired crack in 
   * the k-th iteration is given as crack_low = 0.99^k * 0.3. 
   * We iterate 21 times, thus in the final update we have 
   * crack_low = 0.99^21 * 0.3 = 0.243. 
   * In order to reduce the computational time, here, we use the update
   * crack_low = 0.98^k * 0.3 leading to the final desired crack 
   * crack_low = 0.98^10 * 0.3 = 0.245. 
   * 
  **/
	for (unsigned int i = 0; i < 11; i++)
	{
		try
		{
      // Write output
		  Alg.GetOutputHandler()->Write(out,1);
			LFunc.update_Crack(crack_low_x, crack_high_x, crack_low_y, crack_high_y);
			FUNC4.update_Crack(crack_low_x, crack_high_x, crack_low_y, crack_high_y);

      // Solve the reduced optimal control problem
			Alg.Solve(q);

      // Update desired crack
      crack_low_x *= 0.98;
			out<<"UPDATED DESIRED CRACK Nr. "<<i<< ": varphi_d := ["<<crack_low_x<<","<<
			crack_high_x<<"] x ["<<crack_low_y<<","<<crack_high_y<<"]."<<"\n";
			Alg.GetOutputHandler()->Write(out,1);
      Alg.GetOutputHandler()->ReInit();

		}
    catch (DOpEException &e)
		{
      // Handle exceptions
			std::cout<< "Warning: During execution of `" + e.GetThrowingInstance()
		          + "` the following Problem occurred!" << std::endl;
			std::cout << e.GetErrorMessage() << std::endl;
		}
  }
  return 0;
}
#undef FDC
#undef CDC
#undef FE
#undef DOFHANDLER
