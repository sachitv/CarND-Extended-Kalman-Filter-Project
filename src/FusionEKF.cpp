#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static float const POSITION_UNCERTAINTY = 1;
static float const VELOCITY_UNCERTAINTY = 1000;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  m_PreviousTimestamp = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    ekf_.SetR_Laser(R_laser_);
    ekf_.SetR_Radar(R_radar_);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    m_PreviousTimestamp = measurement_pack.timestamp_;

    MatrixXd P(4, 4);
    P <<    POSITION_UNCERTAINTY, 0, 0, 0,
            0, POSITION_UNCERTAINTY, 0, 0,
            0, 0, VELOCITY_UNCERTAINTY, 0,
            0, 0, 0, VELOCITY_UNCERTAINTY;

    MatrixXd Finit(4, 4);
    Finit <<  1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;

    MatrixXd QInit(4, 4);
    QInit <<  0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0;

    MatrixXd H(2, 4);
    H <<  1, 0, 0, 0,
          0, 1, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      float const rho = measurement_pack.raw_measurements_[0];
      float const phi = measurement_pack.raw_measurements_[1];
      float const rhodot = measurement_pack.raw_measurements_[2];

      float const px = rho * cos(phi);
      float const py = rho * sin(phi);
      float const vx = rhodot * cos(phi);
      float const vy = rhodot * sin(phi);

      VectorXd x(4);
      x << px, py, vx, vy;

      ekf_.Init(x, P, Finit, H, QInit);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      VectorXd x(4);
        float const px = measurement_pack.raw_measurements_(0);
        float const py = measurement_pack.raw_measurements_(1);
      x << px, py, 0, 0;

      ekf_.Init(x, P, Finit, H, QInit);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float const deltaTime = ( measurement_pack.timestamp_ - m_PreviousTimestamp ) / 1000000.0;
  m_PreviousTimestamp = measurement_pack.timestamp_;
  ekf_.Predict(deltaTime);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "m_x = " << ekf_.GetX() << endl;
  cout << "m_P = " << ekf_.GetP() << endl;
}
