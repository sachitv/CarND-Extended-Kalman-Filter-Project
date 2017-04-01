#include "kalman_filter.h"
#include <cassert>

using Eigen::MatrixXd;
using Eigen::VectorXd;

float const KalmanFilter::s_NoiseAx = 9.0;
float const KalmanFilter::s_NoiseAy = 9.0;

static MatrixXd R_laser(2, 2);
static MatrixXd R_radar(3, 3);

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const VectorXd &x_in, const MatrixXd &P_in, const MatrixXd &F_in,
                        const MatrixXd &H_in, const MatrixXd &Q_in) {
  m_x = x_in;
  m_P = P_in;
  m_F = F_in;
  m_H = H_in;
  m_Q = Q_in;

  //measurement covariance matrix - laser
  R_laser <<    0.0225, 0,
                0, 0.0225;

  //measurement covariance matrix - radar
  R_radar <<    0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;
}

void KalmanFilter::Predict(float const deltaTime) {
  /**
  TODO:
    * predict the state
  */
    m_F(0, 2) = deltaTime;
    m_F(1, 3) = deltaTime;


    float const dt2 = deltaTime * deltaTime;
    float const dt3 = dt2 * deltaTime;
    float const dt4 = dt3 * deltaTime;

    m_Q(0, 0) = dt4 / 4 * s_NoiseAx;
    m_Q(0, 2) = dt3 / 2 * s_NoiseAx;
    m_Q(1, 1) = dt4 / 4 * s_NoiseAy;
    m_Q(1, 3) = dt3 / 2 * s_NoiseAy;
    m_Q(2, 0) = dt3 / 2 * s_NoiseAx;
    m_Q(2, 2) = dt2 * s_NoiseAx;
    m_Q(3, 1) = dt3 / 2 * s_NoiseAy;
    m_Q(3, 3) = dt2 * s_NoiseAy;

    //Assuming that u is a zero array here, that's why i'm skipping it
    m_x = m_F * m_x;
    m_P = m_F * m_P * m_F.transpose() + m_Q;


}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    MatrixXd const y = z - m_H * m_x;

    MatrixXd const HT(m_H.transpose());

    MatrixXd const S = m_H * m_P * HT + R_laser;
    MatrixXd const K = m_P * HT * S.inverse();


    m_x = m_x + K * y;

    long const x_size = m_x.size();
    MatrixXd const I(MatrixXd::Identity(x_size, x_size));

    m_P = (I - K * m_H) * m_P;


}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    MatrixXd const hx( GenerateHx() );
    MatrixXd const y = z - hx;

    MatrixXd const Hj( CalculateJacobian() );

    if(Hj.isZero())
    {
        return;
    }

    MatrixXd const HjT(Hj.transpose());
    MatrixXd const S = ( ( Hj * m_P ) * HjT ) + R_radar;
    MatrixXd const SInv(S.inverse());
    MatrixXd const K = ( m_P * HjT ) * SInv;

    m_x = m_x + K * y;

    long const x_size = m_x.size();
    MatrixXd const I(MatrixXd::Identity(x_size, x_size));

    m_P = (I - K * Hj) * m_P;


}

MatrixXd KalmanFilter::CalculateJacobian() const
{
    MatrixXd Hj(3,4);

    Hj <<   0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;

    //recover state parameters
    float const px = m_x(0);
    float const py = m_x(1);
    float const vx = m_x(2);
    float const vy = m_x(3);

    //check division by zero
    if( (px == 0) && (py == 0) )
    {
        //assert( (px!=0) && (py!= 0) );
        return Hj;
    }

    //compute the Jacobian matrix
    float const pxpy = px*px + py*py;

    Hj(0, 0) = px / sqrt(pxpy);
    Hj(0, 1) = py / sqrt(pxpy);
    Hj(0, 2) = 0;
    Hj(0, 3) = 0;

    Hj(1, 0) = - py / (pxpy);
    Hj(1, 1) = px / (pxpy);
    Hj(1, 2) = 0;
    Hj(1, 3) = 0;

    Hj(2, 0) = py * (vx * py - vy * px) / pow(pxpy, 3.0/2.0);
    Hj(2, 1) = px * (vy * px - vx * py) / pow(pxpy, 3.0/2.0);
    Hj(2, 2) = px / sqrt(pxpy);
    Hj(2, 3) = py / sqrt(pxpy);

    return Hj;
}

Eigen::MatrixXd KalmanFilter::GenerateHx() const
{
    MatrixXd hx(3, 1);

    float const px = m_x(0);
    float const py = m_x(1);
    float const vx = m_x(2);
    float const vy = m_x(3);

    float const rho = sqrt((px * px) + (py * py));
    float const phi = atan2(py, px);
    float const rhodot = ((px * vx) + (py * vy)) / rho;

    hx << rho, phi, rhodot;

    return hx;
}

void KalmanFilter::SetR_Radar(const MatrixXd &m_R_Radar) {
    KalmanFilter::m_R_Radar = m_R_Radar;
}

void KalmanFilter::SetR_Laser(const MatrixXd &m_R_Laser) {
    KalmanFilter::m_R_Laser = m_R_Laser;
}
