#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;//false;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;//false;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);
    P_ << MatrixXd::Identity(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = M_PI*8.0/(4.0*4.0); // bike acceleration

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 2*M_PI/(4.0*4.0); // bike turn rate

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

    /**
    TODO:

    Complete the initialization. See ukf.h for other member properties.

    Hint: one or more values initialized above might be wildly off...
    */

    is_initialized_ = false;

    n_x_ = 5;

    n_aug_ = 7;

    lambda_ = 3 - n_aug_;

    weights_ = VectorXd(2*n_aug_ + 1);
    weights_.segment(1, 2*n_aug_).fill(0.5/(n_aug_+lambda_));
    weights_(0) = lambda_ / (lambda_ + n_aug_);

    H_laser_ = MatrixXd(2,5);
    H_laser_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0;

    R_laser_ = MatrixXd(2,2);
    R_laser_ << std_laspx_*std_laspx_,0,
            0,std_laspy_*std_laspy_;

    time_us_ = 0;

    Zsig_ = MatrixXd(3, 2*n_aug_+1);

    z_pred = VectorXd(3);

    S = MatrixXd(3,3);

    R_radar_ = MatrixXd(3,3);
    R_radar_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;

    x_aug = VectorXd(7);
    P_aug = MatrixXd(7, 7);

    Q_ = MatrixXd(2,2);
    Q_ << std_a_ * std_a_, 0,
            0, std_yawdd_ * std_yawdd_;

    Xsig_ = MatrixXd(7, 2*n_aug_+1);
    Xsig_pred_ = MatrixXd(5, 2*n_aug_+1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */

    /*****************************************************************************
     *  Initialization
     ****************************************************************************/

    if (!is_initialized_) {

        // first measurement
        cout << "UKF: " << endl;
        double px = 0.0;
        double py = 0.0;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

            double rho = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];

            px = rho * cos(phi);
            py = rho * sin(phi);

            if(fabs(px) < 0.0001) px = 0.001;
            if(fabs(py) < 0.0001) py = 0.001;

            cout << "(px,py) from RADAR: (" << px << ", " << py << ")"<< endl;

        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

            px = meas_package.raw_measurements_[0];
            py = meas_package.raw_measurements_[1];

            if(fabs(px) < 0.0001) px = 0.001;
            if(fabs(py) < 0.0001) py = 0.001;

            cout << "(px,py) from LIDAR: (" << px << ", " << py << ")" << endl;

            P_(0,0) = 0.15;
            P_(1,1) = 0.15;
        }

        // px, py, v, yaw(psi), yaw(psi)_dot
        x_ << px, py, 0, 0, 0;

        time_us_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }


    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    Prediction(dt);

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (!(x_(0) == 0 && x_(1) == 0)){
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            // Radar update
            if (use_radar_ == true)
                UpdateRadar(meas_package);

        } else {
            // Laser update
            if (use_laser_ == true)
                UpdateLidar(meas_package);
        }
    }

    // print state and covariance
    // cout << "x_ = " << x_ << endl;
    // cout << "P_ = " << P_ << endl;
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

    // augmented x

    x_aug << x_.array(), 0.0, 0.0;

    P_aug.topLeftCorner(5,5) = P_;
    P_aug.bottomRightCorner(2,2) = Q_;

    // square root of P_aug

    MatrixXd A = P_aug.llt().matrixL();

    // sigma points Xsig_

    Xsig_.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++){
        Xsig_.col(i + 1) = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
        Xsig_.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
    }

    // predict sigma points Xsig_pred_

    for (int i = 0; i < 2*n_aug_+1; i++) {
        double px = Xsig_(0, i);
        double py = Xsig_(1, i);
        double v = Xsig_(2, i);
        double yaw = Xsig_(3, i);
        double yawd = Xsig_(4, i);
        double nu_a = Xsig_(5, i);
        double nu_yawdd = Xsig_(6, i);

        if (fabs(yawd) > 0.001) {
            px += v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py += v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        } else {
            px += v * delta_t * cos(yaw);
            py += v * delta_t * sin(yaw);
        }

        v += nu_a * delta_t;
        yaw += yawd * delta_t;

        // add noise
        px += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        yaw += 0.5 * nu_yawdd * delta_t * delta_t;
        yawd += nu_yawdd * delta_t;

        Xsig_pred_.col(i) << px, py, v, yaw, yawd;
    }

    // update new state x_
    x_.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; i++){
        x_ += weights_(i) * Xsig_pred_.col(i);
    }

    // update new state covariance matrix P_
    P_.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; i++){
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        while(x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
        while(x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

        P_ += weights_(i) * x_diff * x_diff.transpose();
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

    VectorXd z = meas_package.raw_measurements_;

    VectorXd y = z - H_laser_ * x_;
    MatrixXd Ht = H_laser_.transpose();
    MatrixXd PHt = P_ * Ht;
    MatrixXd S = H_laser_ * PHt + R_laser_;
    MatrixXd Si = S.inverse();
    MatrixXd K = PHt * Si;
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    // new state and uncertainty covariance.
    x_ = x_ + K * y;
    P_ = (I - K * H_laser_) * P_;

    NIS_laser_ = y.transpose() * Si * y;
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

    // project sigma points (Xsig_pred_) into measurement space (Zsig_)

    for (int i = 0; i < 2*n_aug_ + 1; i++){
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);

        double vx = sin(yaw)*v;
        double vy = cos(yaw)*v;

        double rho = sqrt(px*px+py*py);
        double phi = atan2(py,px);
        double rho_dot = (px*vy+py*vx)/rho;

        Zsig_.col(i) << rho, phi, rho_dot;
    }


    // mean of sigma points Zsig_

    z_pred.fill(0.0);
    for(int i = 0; i < 2*n_aug_+1; i++){
        z_pred += weights_(i) * Zsig_.col(i);
    }

    // covariance of sigma points Zsig_

    S.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; i++){
        VectorXd z_diff = Zsig_.col(i) - z_pred;

        while(z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while(z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    S += R_radar_;


    // cross correlation matrix Tc

    MatrixXd Tc = MatrixXd(n_x_, 3);
    Tc.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; i++){
        // residual
        VectorXd z_diff = Zsig_.col(i) - z_pred;

        while(z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while(z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        while(x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while(x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }


    // Kalman gain K

    MatrixXd K = Tc * S.inverse();


    // residual

    VectorXd y = meas_package.raw_measurements_ - z_pred;

    while(y(1) > M_PI) y(1) -= 2*M_PI;
    while(y(1) < -M_PI) y(1) += 2*M_PI;

    // update state vector and covariance P_

    x_ = x_ + K * y;
    P_ = P_ - K * S * K.transpose();

    NIS_radar_ = y.transpose() * S.inverse() * y;

}
