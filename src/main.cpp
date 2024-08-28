#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "MPC.h"
#include "json.hpp"
#include <fstream>

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

/**
* Create a homogenious transformation matrix that transforms a homogenious 
* 2D point from the vehicle coordinate system to the map coordinate system.
* mapPoint = A . carPoint					
* [mx]   [ cos(psi),	-sin(psi),	tx ]   [cx]
* [my] = [ sin(psi),	 cos(psi),	ty ] . [cy]
* [ 1]   [        0,            0,   1 ]   [ 1]
* ( carPoint = A^(-1) . mapPoint )
* @param psi Rotation of the vehicle CS w.r.t. the map CS
* @param tx Vehicle position along the x-axis of the map CS
* @param ty Vehicle position along the y-axis of the map CS
*/
Eigen::Matrix3d getCar2MapTransformationMatrix(double psi, double tx, double ty) {
  Eigen::Matrix3d A;
  A << cos(psi), -sin(psi), tx, sin(psi), cos(psi), ty, 0., 0., 1.; 
  return A;
}

/**
* Test function for the transformation matrix.
*/
void myTestFunction() {
  // Test coordinate transformation
  std::cout << "Testing coordinate transformations:" << std::endl;
  Eigen::Matrix3d A = getCar2MapTransformationMatrix(pi()/4., 1., 3.);
  std::cout << "Car to map CS:\n" << A << std::endl;
  std::cout << "Map to car CS:\n" << A.inverse() << std::endl;
  
  // Transform point from car CS to map CS and back
  Eigen::Vector3d point_car(1., 0., 1.);
  Eigen::Vector3d point_map = A * point_car;
  std::cout << "Point in car CS:\n" << point_car << std::endl;
  std::cout << "Point in map CS: mp=A.cp\n" << point_map << std::endl;
  std::cout << "Back to car CS: cp=A^(-1).mp\n" << A.inverse() * point_map << std::endl;
}


/** 
* Calculate polynomial coefficients represented in the vehicle 
* coordinate system. 
* @param tx Vehicle position along the x-axis of the map CS
* @param ty Vehicle position along the y-axis of the map CS
* @param psi Rotation of the vehicle CS w.r.t. the map CS
* @param ptsx Vector of x-values of the points to be fitted 
* @param ptsx Vector of y-values of the points to be fitted
* @param poly_order Order of the polynomial
*/
Eigen::VectorXd getPolyCoeffs(double px, double py, double psi, std::vector<double> ptsx, std::vector<double> ptsy, int poly_order) {
  // transform waypoints to car coordinate system
  Eigen::Matrix3d A_c2m = getCar2MapTransformationMatrix(psi, px, py);	// A_c2m transforms from car CS to map CS
  Eigen::Matrix3d A_m2c = A_c2m.inverse();	// A_m2c transforms from map CS to car CS
  
  std::vector<double> ptsx_car;
  std::vector<double> ptsy_car;
  for (int i=0; i<ptsx.size(); i++) {
    Eigen::Vector3d point_map(ptsx[i], ptsy[i], 1.);
    Eigen::Vector3d point_car = A_m2c * point_map;
    ptsx_car.push_back(point_car[0]);
    ptsy_car.push_back(point_car[1]);
  }
  
  // copy std::vector data to Eigen::VectorXd
  Eigen::Map<Eigen::VectorXd> eigen_ptsx(ptsx_car.data(), ptsx_car.size());
  Eigen::Map<Eigen::VectorXd> eigen_ptsy(ptsy_car.data(), ptsy_car.size());
  
  // fit polynomial to the transformed waypoints
  Eigen::VectorXd poly_coeffs = polyfit(eigen_ptsx, eigen_ptsy, poly_order);
  return poly_coeffs;
}
  
/** 
* Calculate the state vector based on the telemetry provided by the simulation.
* @param x Position of the vehicle along the x-axis of the map
* @param y Position of the vehicle along the y-axis of the map
* @param psi Rotation of the vehicle CS w.r.t. the map CS
* @param v Velocity of the vehicle
* @param delta Steering angle of the vehicle
* @param a Acceleration of the vehicle
* @param poly_coeffs Coefficients of the polynomial representing the desired vehicle trajectory 
*/
Eigen::VectorXd getStateVector(double x, double y, double psi, double v, double delta, double a, Eigen::VectorXd  poly_coeffs) {
  // calculate the cross track error 
  double cte = polyeval(poly_coeffs, 0.0);

  // calculate the orientational error
  double epsi = -atan(poly_coeffs[1]);
  
  // coordinates w.r.t. the car CS
  double x_car = 0;
  double y_car = 0;
  // rotation of the vehicle w.r.t. the target trajectory polynomial 
  // at the origin of the vehicle CS is zero
  double psi_car = 0; 
  
  Eigen::VectorXd state(8);
  // steering angle delta has the opposite sign in the simulation 
  state << x_car, y_car, psi_car, v, cte, epsi, -delta, a;
  return state;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;


  // Read parameters file
  std::ifstream config_file("../src/parameters.json");
  if (!config_file.is_open()) {
    std::cerr << "Error: Could not open parameters file!" << std::endl;
    return 1;
  } else {
    std::cout << "Parameters file opened successfully!" << std::endl;
  }

  nlohmann::json config;
  config_file >> config;

  double cost_factor_cte = config["cost_factor_cte"];
  double cost_factor_epsi = config["cost_factor_epsi"];
  double cost_factor_v = config["cost_factor_v"];
  double v_ref = config["v_ref"];

  cout << "cost_factor_cte: " << cost_factor_cte << endl;
  cout << "cost_factor_epsi: " << cost_factor_epsi << endl;
  cout << "cost_factor_v: " << cost_factor_v << endl;
  cout << "v_ref: " << v_ref << endl;

  mpc.SetCostFactors(cost_factor_cte, cost_factor_epsi, cost_factor_v);
  mpc.SetVRef(v_ref);
  
  //myTestFunction();

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    // cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          // adding actuator state parameters
          double delta = j[1]["steering_angle"];
		  double a = j[1]["throttle"];
         
          // Set up transformation matrix 
          Eigen::Matrix3d A_c2m = getCar2MapTransformationMatrix(psi, px, py);	// transforms from car CS to map CS
          Eigen::Matrix3d A_m2c = A_c2m.inverse();	// transforms from map CS to car CS
         
          /*
          * Calculate steering angle and throttle using MPC.
          * Both are in between [-1, 1].
          */
          // Fit polynomial to the waypoints
          int poly_order = 3;
          auto poly_coeffs = getPolyCoeffs(px, py, psi, ptsx, ptsy, poly_order);
          
          // Calculate the state vector
          auto state = getStateVector(px, py, psi, v, delta, a, poly_coeffs);

          // Call mpc
          auto vars = mpc.Solve(state, poly_coeffs);
     
          double steering_angle;
          double throttle_value;
          
          // steering angle has opposite sign in the simulation
          steering_angle = (-vars[vars.size()-2] / deg2rad(25));
          throttle_value = vars[vars.size()-1]; 
         

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steering_angle;
          msgJson["throttle"] = throttle_value;

          // Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          // add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line
          
          // The vars vector returned by the mpc contains the x- and y-values
          // of the optimized vehicle trajectory and the actuator parameters
          size_t n_points = (vars.size()-2) / 2;
          //std::cout << "n_points: " << n_points << std::endl;
          for (int i=1; i<n_points; i++) {
            mpc_x_vals.push_back(vars[i]);
            mpc_y_vals.push_back(vars[i+n_points]);
          }

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          // Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          
          // add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
          
          // Transform waypoints to the vehicle coordinate system
          
          // Generate waypoints from the polynomial
          Eigen::Vector3d point_max_dist_map(ptsx[ptsx.size()-1], ptsy[ptsy.size()-1], 1.);
          Eigen::Vector3d point_max_dist_car = A_m2c * point_max_dist_map;

          const size_t waypoint_cnt = 20;
          for (int i=1; i<waypoint_cnt; i++) {
            double x_val = point_max_dist_car[0] / waypoint_cnt * i;
            double y_val = polyeval(poly_coeffs, x_val);
            next_x_vals.push_back(x_val);
            next_y_vals.push_back(y_val); 
          }
 
          
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
