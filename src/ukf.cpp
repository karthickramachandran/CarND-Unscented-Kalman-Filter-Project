#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3.0;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 2.0;

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

    // No of elements, state dimension
    n_x_ = 5;

    // No of elements inkl noise, augmented dimension
    n_aug_ = 7;

    //time when the state is true, in us
    time_us_ = 0;

    // initialization flag
    is_initialized_ = false;

    //predicted sigma points matrix
    Xsig_pred_ = Eigen::MatrixXd(n_x_, 2*n_aug_+1);

    //weights w
    weights_ = Eigen::VectorXd(2*n_aug_+1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    //1. Part-1: If not initialized
    //            -> Init LiDAR
    //            -> Init Radar
    //2. Part 2: if already initializated

    if(!is_initialized_)
    {
        // Initial state vector
        x_<<1,1,1,1,1;

        // initial covariance matrix 5x5
        P_ <<   1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1 ;

        // Lidar has px, py and the rest v, yaw, yawd could be set to zero
        if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            x_(0)= meas_package.raw_measurements_(0);
            x_(1)= meas_package.raw_measurements_(1);
            x_(2)= 0;
            x_(3)= 0;
            x_(4)= 0;
        }

        else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            // The radar measure ment contains only three elements
            // rho, phi and rho_dot

            double rho      = meas_package.raw_measurements_(0);
            double phi      = meas_package.raw_measurements_(1);
            double rho_dot  = meas_package.raw_measurements_(2);

            // Convert from polar to cartesian
            // px=rho*cos(phi)
            // py=rho*sin(phi)
            x_(0) = rho*cos(phi);
            x_(1) = rho*sin(phi);
            //
            x_(0) = sqrt(rho_dot*cos(phi)+rho_dot*sin(phi));
            x_(0) = 0;
            x_(0) = 0;
        }

        // Set the is_initialized_ flag to true
        is_initialized_ = true;

        // Copy the timestamp
        time_us_ = meas_package.timestamp_;

        return;
    }

    // If already initialized...
    // Update teh measurement depending upon the sensor type(LiDAR or RADAR)

    double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_ ;

    Prediction(delta_t);

    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        UpdateLidar(meas_package);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        UpdateRadar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

    /*Step1: Generating Sigma points*/
    // Sigma Point Matrix

    lambda_ = 3 - n_x_;
    Eigen::MatrixXd Xsig = Eigen::MatrixXd(n_x_, 2*n_aug_+1);

    // Square root of P
    Eigen::MatrixXd A = P_.llt().matrixL();

    Xsig.col(0) = x_;

    for (int i = 0; i < n_x_; i++)
    {
        Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
        Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
    }

    /*Augmentation*/
    lambda_ = 3 - n_aug_;
    //Augmented Mean vector
    Eigen::VectorXd x_aug = Eigen::VectorXd(n_aug_); //n_aug_ is 7

    //create augmented state covariance

    Eigen::MatrixXd P_aug = Eigen::MatrixXd(n_aug_, n_aug_); //7X7
    // Augmented sigma point matrix
    Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd(n_aug_, 2 * n_aug_ + 1);

    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    //create square root matrix
    Eigen::MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }

    /*Predict*/
    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001)
        {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else
        {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }

    //create vector for weights
    //Eigen::VectorXd weights = Eigen::VectorXd(2*n_aug_ + 1);

    // set weights
    double weight_0 = lambda_/(lambda_ + n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) //2n+1 weights
    {
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }

    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) //iterate over sigma points
    {
        x_ = x_+ weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) //iterate over sigma points
    {
        // state difference
        Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    Eigen::VectorXd z = meas_package.raw_measurements_;

    //set measurement dimension, LiDAR can measure Px and Py
    int n_z = 2;


    Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * n_aug_ + 1);

    for (int i=0; i< 2 * n_aug_ + 1; i++)
    {
        Zsig(0, i) = Xsig_pred_(0,i);
        Zsig(1, i) = Xsig_pred_(1,i);
    }

    Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
    z_pred.fill(0.0);

    for (int i=0; i < 2*n_aug_+1; i++)
    {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    Eigen::MatrixXd S = Eigen::MatrixXd(n_z,n_z);

    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        //residual
        Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    Eigen::MatrixXd R = Eigen::MatrixXd(n_z,n_z);
    R <<    std_laspx_*std_laspx_   , 0,
            0                       , std_laspy_*std_laspy_;
    S = S + R;


    //calculate cross correlation matrix
    Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++)   //2n+1 simga points
    {

        //residual
        Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

    Eigen::VectorXd z = meas_package.raw_measurements_;

    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //create matrix for sigma points in measurement space
    Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        double px = Xsig_pred_(0, i);
        double py = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double psi = Xsig_pred_(3, i);

        double rho = sqrt(px*px+py*py);
        double pi= atan2(py, px);
        double rho_dot=(px*cos(psi)*v+py*sin(psi)*v)/rho;

        Zsig(0, i)=rho;
        Zsig(1, i)=pi;
        Zsig(2, i)=rho_dot;
    }

    //mean predicted measurement
    Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
    z_pred.fill(0.0);

    for (int i=0; i < 2*n_aug_+1; i++)
    {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //innovation covariance matrix S
    Eigen::MatrixXd S = Eigen::MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        //residual
        Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    Eigen::MatrixXd R = Eigen::MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_,    0,                          0,
            0,                      std_radphi_*std_radphi_,    0,
            0,                      0,                          std_radrd_*std_radrd_;
    S = S + R;

    //create matrix for cross correlation Tc
    Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        //residual
        Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;

        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    Eigen::MatrixXd K = Tc * S.inverse();

    //residual
    Eigen::VectorXd z_diff = z - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    //update state mean and covariance matrix
    x_= x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
}
