#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <limits>


// #include "param.h"

using CppAD::AD;

// Set the timestep length and duration
size_t N = 25;
double dt = 0.05;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// NOTE: feel free to play around with this
// or do something completely different
// double v_ref = 40; //65; 

// The solver takes all the state variables and actuator
// variables in a single vector. Thus, we should establish
// when one variable starts and another ends to make our lifes easier.
const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { 
    this->coeffs = coeffs; 
    cost_factor_cte = 1;
    cost_factor_epsi = 1;
    cost_factor_v = 1;  
  }

  void SetCostFactors(double cost_factor_cte, double cost_factor_epsi, double cost_factor_v)
  {
    this->cost_factor_cte = cost_factor_cte;
    this->cost_factor_epsi = cost_factor_epsi;
    this->cost_factor_v = cost_factor_v;
  }
  void SetVRef(double v_ref)
  {
    this->v_ref = v_ref;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // Implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    
    // The cost is stored in the first element of `fg`. 
    fg[0] = 0;
    
    // Reference State Cost
    // adding cross tracking and orientation error
    for (int t=0; t<N; t++) {
   	  fg[0] += cost_factor_cte * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += cost_factor_epsi * CppAD::pow(vars[epsi_start + t], 2);
      // fg[0] += cost_factor_v * CppAD::pow(vars[v_start + t] - v_ref, 2);
      fg[0] += cost_factor_v * CppAD::pow(vars[v_start + t] - v_ref, 2);
    }
    
    // adding actuator to cost function
    for (int t=0; t<N-1; t++) {
   	  fg[0] += 50 * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += 5 * CppAD::pow(vars[a_start + t], 2);
    }
    
    // minimizing value gap between sequential actuations
    for (int t=0; t<N-2; t++) {
   	  fg[0] += 500 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += 50 * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }
    
    // Setup Constraints

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];
    
    // The rest of the constraints
    for (int t = 1; t < N; t++) {
	  // state parameters      
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];
      
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // control parameters     
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];
      
      // AD<double> f0 = coeffs[0] + coeffs[1] * x0;
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      //AD<double> psides0 = CppAD::atan(coeffs[1]);
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0,2));
      
      // Setup the rest of the model constraints
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt); 
      fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + (v0 * delta0 / Lf * dt));
    } 
  }

private:
  double cost_factor_cte;
  double cost_factor_epsi;
  double cost_factor_v;
  double v_ref;

};

//
// MPC class definition implementation.
//
MPC::MPC() :
  cost_factor_cte(1),
  cost_factor_epsi(1),
  cost_factor_v(1),
  v_ref(1) 
  {}

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];
  
  // actuator state parameters
  double delta = state[6];
  double a = state[7];

  // Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  //
  // N 	... number of timesteps
  // dt	... timestep
  // T	... time horizon (T = N * dt)
  // state 		... [x, y, psi, cte, epsi] (vector of size 6)
  // actuation	... [delta, a] (vector of size 2)
  size_t n_state_params = 6;
  size_t n_actuation_params = 2;  
  size_t n_vars = n_state_params * N + (N - 1) * n_actuation_params;
  //std::cout << "n_vars: " << n_vars << std::endl;
  
  // Set the number of constraints
  //
  // we want to constrain the state parameters since these are governed by the model dynamics
  // (the actuation parameters are to be optimized and therefore not constrained)
  size_t n_constraints = n_state_params * N;
  
  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  
  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi; 
  
  // Set the initial values of the actuators
  vars[delta_start] = delta;
  vars[a_start] = a;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Set lower and upper limits for variables.
  
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = std::numeric_limits<double>::lowest();
    vars_upperbound[i] = std::numeric_limits<double>::max();
  }
  
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  // Constrain the steering actuation of the first two time step to model the latency
  vars_lowerbound[delta_start] = delta;
  vars_upperbound[delta_start] = delta;
  vars_lowerbound[delta_start + 1] = delta;
  vars_upperbound[delta_start + 1] = delta;
  
  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  // Constrain the throttle actuation of the first two time step to model the latency
  vars_lowerbound[a_start] = a;
  vars_upperbound[a_start] = a;
  vars_lowerbound[a_start + 1] = a;
  vars_upperbound[a_start + 1] = a;
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;
  
 
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);
  fg_eval.SetCostFactors(cost_factor_cte, cost_factor_epsi, cost_factor_v);
  fg_eval.SetVRef(v_ref);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  //std::vector<double> x_vals(solution.x
  
  // Copy x- and y-values of the optimal trajectory
  const size_t n_elems = N*2 + 2;
  std::vector<double> output(n_elems);
  for (int i=0; i<N; i++) {
    output[i]   = solution.x[x_start + i];
    output[i+N] = solution.x[y_start + i];
  }
  
  // Return the unconstrained actuation parameters at time step (t0 + 2*dt) 
  // which equals (t0 + 100 ms) 
  output[2*N]   = solution.x[delta_start+2];
  output[2*N+1] = solution.x[a_start+2];
  
  return output; 
}

void MPC::SetCostFactors(double cost_factor_cte, double cost_factor_epsi, double cost_factor_v)
{
  this->cost_factor_cte = cost_factor_cte;
  this->cost_factor_epsi = cost_factor_epsi;
  this->cost_factor_v = cost_factor_v;
}

void MPC::SetVRef(double v_ref)
{
  this->v_ref = v_ref;
}
