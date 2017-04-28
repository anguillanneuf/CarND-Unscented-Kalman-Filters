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
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

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

    H_laser_ = MatrixXd(2,5);
    H_laser_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0;

    R_laser_ = MatrixXd(2,2);
    R_laser_ << std_laspx_*std_laspx_,0,
        0,std_laspy_*std_laspx_;

    time_us_ = 0.0;

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
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "UKF: " << endl;
        float px = 0.0;
        float py = 0.0;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Extract px, py from radar's polar coordinates.
            */
            float rho = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];

            px = rho * cos(phi);
            py = rho * sin(phi);

            cout << "(px,py) from RADAR: " << px << ", " << py << endl;

        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {


            px = meas_package.raw_measurements_[0];
            py = meas_package.raw_measurements_[1];

            // is this step necessary?
            if(fabs(px==0.0)) px=1.0e-3;
            if(fabs(py==0.0)) py=1.0e-3;

            cout << "(px,py) from LIDAR: " << px << ", " << py << endl;
        }

        time_us_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }


    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    Prediction(dt);


    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar update
        if(use_radar_== true)
            UpdateRadar(meas_package);

    } else {
        // Laser update
        if (use_laser_ == true)
            UpdateLidar(meas_package);
    }

    // print the output
    cout << "x_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;
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
}
