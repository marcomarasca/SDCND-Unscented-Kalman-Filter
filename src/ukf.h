#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "tools.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <functional>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  // Initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // if this is false uses the UKF equations for the update step, otherwise uses
  // the standard KF equations (more efficient)
  bool use_laser_kf;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  // state covariance matrix
  MatrixXd P_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;

  // Latest NIS (Normalized Innovation Squared) computed for the radar
  double NIS_radar_;

  // Latest NIS (Normalized Innovation Squared) computed for the laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage &measurement_pack);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage &measurement_pack);

private:

  // Radar measurements dimension: rho, phi, rho_dot
  int _n_z_radar = 3;

  // Laser measurements dimentsion: px, py
  int _n_z_laser = 2;

  // predicted sigma points matrix
  MatrixXd _Xsig_pred;

  // measurement covariance matrix - radar
  MatrixXd _R_radar;

  // measurement covariance matrix - laser
  MatrixXd _R_laser;

  // measurement matrix for laser
  MatrixXd _H_laser;

  // Identity matrix
  MatrixXd _I;

  // Weights of sigma points
  VectorXd _weights;

  // Tool object used to work on polar coordinates
  Tools _tools;

  /**
   * Initilizes the kalman filter using the given measurement package, if the data comes from a radar
   * sensor converts the state coordinates from polar to cartesian
   * 
   * @param meas_package The measurement at k+1
   */
  void _InitializeFilter(const MeasurementPackage &measurement_pack);

  /**
   * Generates sigma points for the UKF
   * 
   * @return The matrix with the generated (augmented) sigma points
   */
  MatrixXd _GenerateSigmaPoints();

  /**
   * Updates the sigma points prediction for the generated sigma points, given the time span
   * 
   * @param Xsig_aug The matrix of generated sigma points
   */
  void _UpdateSigmaPointsPrediction(const MatrixXd &Xsig_aug, double delta_t);

  /**
   * Updates the mean state and covariance given the sigma points prediction
   */
  void _UpdateMeanAndCovariancePrediction();

  /**
   * Updates the state and the state covariance matrix given the error, measurement covariance matrix and
   * the cross correlation matrix, updates the given NIS value
   * 
   * @param z_diff Difference between sensor measurement and prediction
   * @param S Measurement covariance matrix
   * @param Tc Matrix for cross correlation
   * @param NIS_out NIS out
   */ 
  void _UpdateStateUKF(const VectorXd &z_diff, const MatrixXd &S, const MatrixXd &Tc, double &NIS_out);

  /**
   * Function pointer for the lidar update that is set to either _UpdateLidarKF or _UpdateLidarUKF according
   * to the use_laser_kf value.
   */
  void (UKF::*_UpdateLidar)(const MeasurementPackage &measurement_pack);

  /**
   * Updates the mean state and covariance matrix using UKF (sigma points) 
   * 
   * @param measurement_package The measurement at k+1
   */
  void _UpdateLidarUKF(const MeasurementPackage &measurement_pack);

  /**
   * Updates the mean state and covariance matrix using standard KF
   * 
   * @param measurement_package The measurement at k+1
   */
  void _UpdateLidarKF(const MeasurementPackage &measurement_pack);

};

#endif /* UKF_H */
