#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  void SetCostFactors(double cost_factor_cte, double cost_factor_epsi, double cost_factor_v);
  void SetVRef(double v_ref);

  double getCostFactorCte() const { return cost_factor_cte; }
  double getCostFactorEpsi() const { return cost_factor_epsi; }
  double getCostFactorV() const { return cost_factor_v; }
  double getVRef() const { return v_ref; }  

private:
  double cost_factor_cte;
  double cost_factor_epsi;
  double cost_factor_v;
  double v_ref;

};

#endif /* MPC_H */
