#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:
  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param H_in Measurement matrix
   * @param R_in Measurement covariance matrix
   * @param Q_in Process covariance matrix
   */
  void Init(const Eigen::VectorXd &x_in, const Eigen::MatrixXd &P_in, const Eigen::MatrixXd &F_in,
            const Eigen::MatrixXd &H_in, const Eigen::MatrixXd &Q_in);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict(float const deltaTime);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

  void SetR_Radar(const Eigen::MatrixXd &m_R_Radar);
  void SetR_Laser(const Eigen::MatrixXd &m_R_Laser);

  Eigen::MatrixXd GetX() const { return m_x; }
  Eigen::MatrixXd GetP() const { return m_P; }

  static float const s_NoiseAx;
  static float const s_NoiseAy;

private:

  Eigen::MatrixXd GenerateHx() const;
  Eigen::MatrixXd CalculateJacobian() const;

  // state vector
  Eigen::VectorXd m_x;

  // state covariance matrix
  Eigen::MatrixXd m_P;

  // state transistion matrix
  Eigen::MatrixXd m_F;

  // process covariance matrix
  Eigen::MatrixXd m_Q;

  // measurement matrix
  Eigen::MatrixXd m_H;

  Eigen::MatrixXd m_R_Radar;
  Eigen::MatrixXd m_R_Laser;

};

#endif /* KALMAN_FILTER_H_ */
