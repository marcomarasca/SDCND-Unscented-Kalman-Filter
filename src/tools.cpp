#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);

  rmse << 0, 0, 0, 0;

  if (estimations.size() == 0) {
    cout << "Error - The estimations vector cannot be empty.";
    return rmse;
  }

  if (estimations.size() != ground_truth.size()) {
    cout << "Error - The estimations and ground truth vectors should be of the "
            "same size.";
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse /= estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

float Tools::NormalizeAngle(float angle) {
  return atan2(sin(angle), cos(angle));
}