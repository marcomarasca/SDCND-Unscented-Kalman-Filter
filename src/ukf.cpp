#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 0.0001

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  // Initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // If this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // If this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // State dimension
  n_x_ = 5;
  
  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Augmented state dimension
  n_aug_ = 7;

  // Last computed NIS (Normalized Innovation Squared) for radar
  NIS_radar_ = 0.0;

  // Last computed NIS (Normalized Innovation Squared) for radar
  NIS_laser_ = 0.0;

  // Last recorded measurement timestamp
  time_us_ = 0.0;

  // Initial state vector
  x_ = VectorXd::Zero(n_x_);

  // Initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Radar measurement dimension: rho, phi, rho_dot
  _n_z_radar = 3;

  // Laser measurement dimension: px, py
  _n_z_laser = 2;

  // Predicted sigma points
  _Xsig_pred = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  // Weights of sigma points
  _weights = VectorXd::Zero(2 * n_aug_ + 1);

  // Radar measurement covariance matrix
  _R_radar = MatrixXd::Zero(_n_z_radar, _n_z_radar);

  // Laser measurement covariance matrix
  _R_laser = MatrixXd::Zero(_n_z_laser, _n_z_laser);

  // Laser measurement matrix
  _H_laser = MatrixXd::Zero(_n_z_laser, n_x_);

  // Identity Matrix
  _I = MatrixXd::Identity(n_x_, n_x_);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // First measurement, filter initialization
    _InitializeFilter(measurement_pack);
    return;
  }

  // delta time - expressed in seconds
  float delta_t = (measurement_pack.timestamp_ - time_us_) / 1000000.0;

  time_us_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (use_radar_ && measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(measurement_pack);
  } else if (use_laser_ && measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(measurement_pack);
  }

  //cout << "x_ = " << endl << x_ << endl;
  //cout << "P_ = " << endl << P_ << endl;
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // Generate the sigma points
  MatrixXd Xsig_aug = _GenerateSigmaPoints();
  // Predicts the sigma points and updates the prediction
  _UpdateSigmaPointsPrediction(Xsig_aug, delta_t);
  // Updates the state mean and covariance prediction
  _UpdateMeanAndCovariancePrediction();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Update Prediction into Measurement space
   ****************************************************************************/

   // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(_n_z_radar, 2 * n_aug_ + 1);
 
  // Transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px  = _Xsig_pred(0, i);
    double py  = _Xsig_pred(1, i);
    double v   = _Xsig_pred(2, i);
    double yaw = _Xsig_pred(3, i);
    
    double rho     = sqrt(px * px + py * py);
    double phi     = atan2(py, px);
    double rho_dot = 0.0;
    
    if (fabs(rho) > EPS) {
      rho_dot = (px * cos(yaw) * v + py * sin(yaw) * v) / rho;
    }
    
    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    Zsig(2, i) = rho_dot;
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(_n_z_radar);
  
  // Calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred += _weights(i) * Zsig.col(i);
  }

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(_n_z_radar, _n_z_radar);

  // Calculate innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = _tools.NormalizeAngle(z_diff(1));
    S += _weights(i) * z_diff * z_diff.transpose();
  }
  
  S += _R_radar;

  /*****************************************************************************
   *  Update State
   ****************************************************************************/

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, _n_z_radar);

  // Calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    
    VectorXd x_diff = _Xsig_pred.col(i) - x_;
    x_diff(3) = _tools.NormalizeAngle(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = _tools.NormalizeAngle(z_diff(1));
    
    Tc += _weights(i) * x_diff * z_diff.transpose();
  }

  VectorXd z = measurement_pack.raw_measurements_;
  
  VectorXd z_diff = z - z_pred;
  
  z_diff(1) = _tools.NormalizeAngle(z_diff(1));

  _UpdateStateUKF(z_diff, S, Tc, NIS_radar_);
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Update Prediction into Measurement space
   ****************************************************************************/

   // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(_n_z_laser, 2 * n_aug_ + 1);

  // Transform sigma points into measurement space
  Zsig.row(0) = _Xsig_pred.row(0);
  Zsig.row(1) = _Xsig_pred.row(1);
  
  // Mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(_n_z_laser);
  
  // Calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred += _weights(i) * Zsig.col(i);
  }

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(_n_z_laser, _n_z_laser);

  // Calculate innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += _weights(i) * z_diff * z_diff.transpose();
  }
  
  S += _R_laser;

  /*****************************************************************************
   *  Update State
   ****************************************************************************/

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, _n_z_laser);

  // Calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    VectorXd x_diff = _Xsig_pred.col(i) - x_;
    x_diff(3) = _tools.NormalizeAngle(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;

    Tc += _weights(i) * x_diff * z_diff.transpose();
  }

  VectorXd z = measurement_pack.raw_measurements_;
  
  VectorXd z_diff = z - z_pred;

  _UpdateStateUKF(z_diff, S, Tc, NIS_laser_);
}

void UKF::_InitializeFilter(const MeasurementPackage &measurement_pack) {
  string sensor_name = measurement_pack.sensor_type_ == 0 ? "LASER" : "RADAR";

  cout << "Filter Initialization using " << sensor_name << " sensor:" << endl;
  
  float px, py;
  double xy_cov;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar outputs polar coordinates
    float rho = measurement_pack.raw_measurements_(0);
    float phi = measurement_pack.raw_measurements_(1);
    // Converts from polar to cartesian coordinates
    px = rho * cos(phi);
    py = rho * sin(phi);
    xy_cov = std_radr_ * std_radr_;
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Laser outputs raw px and py directly, not conversion necessary
    px = measurement_pack.raw_measurements_(0);  // px
    py = measurement_pack.raw_measurements_(1);  // py
    xy_cov = std_laspx_ * std_laspx_;
  }

  // px, py, v, yaw, yaw_d
  x_ << px, py, 5, 0, 0;

  // Initializes covariance matrix
  P_.setIdentity();
  
  // Set initiale covariance for x and y
  P_(0, 0) = xy_cov;
  P_(1, 1) = xy_cov;
  
  // Initializes the radar covariance matrix
  _R_radar(0, 0) = std_radr_ * std_radr_;
  _R_radar(1, 1) = std_radphi_ * std_radphi_;
  _R_radar(2, 2) = std_radrd_ * std_radrd_;

  // Initializes the laser covariance matrix
  _R_laser(0, 0) = std_laspx_ * std_laspx_;
  _R_laser(1, 1) = std_laspy_ * std_laspy_;
  
  // Initialized the laser measurement matrix
  _H_laser << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;
  
  // set weights
  _weights.fill(0.5 / (lambda_ + n_aug_));
  _weights(0) = lambda_ / (lambda_ + n_aug_);
  
  // Records current timestamp
  time_us_ = measurement_pack.timestamp_;

  // Done initializing
  is_initialized_ = true;

}

MatrixXd UKF::_GenerateSigmaPoints() {
  // Create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  
  // Augmented mean state
  x_aug.head(n_x_) = x_;

  MatrixXd Q = MatrixXd::Zero(n_aug_ - n_x_, n_aug_ - n_x_);
  Q(0, 0) = std_a_ * std_a_;
  Q(1, 1) = std_yawdd_ * std_yawdd_;
  
  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  // A covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q;

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  // Create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  // Create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  
  double lambda_n_aug_sq = sqrt(lambda_ + n_aug_);
  for (int i = 0; i < n_aug_; i++) {
    MatrixXd c1 = lambda_n_aug_sq * L.col(i);
    Xsig_aug.col(i + 1) = x_aug + c1;
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - c1;
  }

  return Xsig_aug;

}

void UKF::_UpdateSigmaPointsPrediction(const MatrixXd &Xsig_aug, double delta_t) {
  // Predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      double px       = Xsig_aug(0, i);
      double py       = Xsig_aug(1, i);
      double v        = Xsig_aug(2, i);
      double yaw      = Xsig_aug(3, i);
      double yaw_d    = Xsig_aug(4, i);
      double nu_a     = Xsig_aug(5, i);
      double nu_yaw_d = Xsig_aug(6, i);
      
      double px_pred, py_pred;
      double sin_yaw = sin(yaw);
      double cos_yaw = cos(yaw);

      // Avoid division by zero
      if (fabs(yaw_d) > EPS) {
        double c1 = v / yaw_d;
        double c2 = yaw + yaw_d * delta_t;
        px_pred = px + c1 * (sin(c2) - sin_yaw);
        py_pred = py + c1 * (cos_yaw - cos(c2));
      } else {
        px_pred = px + v * cos_yaw * delta_t;
        py_pred = py + v * sin_yaw * delta_t;
      }
      
      double v_pred = v; // + 0
      double yaw_pred = yaw + yaw_d * delta_t;
      double yaw_d_pred = yaw_d; // + 0
      
      // Add noise
      double delta_t_2 = 0.5 * delta_t * delta_t;
      
      px_pred    += delta_t_2 * cos_yaw * nu_a;
      py_pred    += delta_t_2 * sin_yaw * nu_a;
      v_pred     += delta_t * nu_a;
      yaw_pred   += delta_t_2 * nu_yaw_d;
      yaw_d_pred += delta_t * nu_yaw_d;
      
      // Updates predicted sigma points into right column
      _Xsig_pred(0, i) = px_pred;
      _Xsig_pred(1, i) = py_pred;
      _Xsig_pred(2, i) = v_pred;
      _Xsig_pred(3, i) = yaw_pred;
      _Xsig_pred(4, i) = yaw_d_pred;
  }

}

void UKF::_UpdateMeanAndCovariancePrediction() {
  
  // Predict state mean
  x_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    x_ += _weights(i) * _Xsig_pred.col(i);
  }
  
  // Predict state covariance matrix
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = _Xsig_pred.col(i) - x_;
    // Normalize angle
    x_diff(3) = _tools.NormalizeAngle(x_diff(3));
    P_ += _weights(i) * x_diff * x_diff.transpose();
  }
}

void UKF::_UpdateStateUKF(const VectorXd &z_diff, const MatrixXd &S, const MatrixXd &Tc, double &NIS_out) {
  MatrixXd Si = S.inverse();

  // Calculate Kalman gain K;
  MatrixXd K = Tc * Si;
  
  // Update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // Computes NIS
  NIS_out = z_diff.transpose() * Si * z_diff;
}
