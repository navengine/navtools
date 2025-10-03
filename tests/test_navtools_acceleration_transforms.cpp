#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <navtools/core/acceleration-transforms.hpp>
#include <navtools/core/coordinate-dcms.hpp>
#include <navtools/core/position-transforms.hpp>
#include <navtools/core/velocity-transforms.hpp>

// Test fixture to hold common test data and pre-calculated values
class AccelerationTransformTest : public ::testing::Test {
 protected:
  // Tolerance for floating-point comparisons
  const double eps_ = 1e-9;

  // Time elapsed for ECI transforms
  const double dt_ = 3600.0;  // 1 hour

  // Reference point (Los Angeles, CA)
  const Eigen::Vector3d lla_la_ = {0.594314087066431, -2.06376826264174, 70.745};
  const Eigen::Vector3d lla0_ = {0.594, -2.063, 0.0};

  // A vector for position, velocity, and acceleration to test with (m, m/s, m/s^2)
  Eigen::Vector3d r_ecef_expected_ = nt::lla2ecef(lla_la_);
  Eigen::Vector3d v_ecef_expected_ = {1.23, -4.56, 7.89};
  Eigen::Vector3d a_ecef_expected_ = {0.123, -0.456, 0.789};

  // Pre-calculated values for comparison (based on the CORRECT formulas)
  Eigen::Vector3d a_ned_expected_ = nt::ecef2nedDcm(lla0_) * a_ecef_expected_;
  Eigen::Vector3d a_enu_expected_ = nt::ecef2enuDcm(lla0_) * a_ecef_expected_;
  Eigen::Vector3d a_eci_expected_ =
      nt::ecef2eciDcm(dt_) *
      (a_ecef_expected_ + 2.0 * nt::OMEGA_ECEF<double>.cross(v_ecef_expected_) +
       nt::OMEGA_ECEF<double>.cross(nt::OMEGA_ECEF<double>.cross(r_ecef_expected_)));
  ;
};

// Test ECEF to NED conversion
TEST_F(AccelerationTransformTest, EcefToNedAcceleration) {
  Eigen::Vector3d a_ned_result = nt::ecef2neda(a_ecef_expected_, lla0_);
  EXPECT_NEAR(a_ned_result(0), a_ned_expected_(0), eps_);
  EXPECT_NEAR(a_ned_result(1), a_ned_expected_(1), eps_);
  EXPECT_NEAR(a_ned_result(2), a_ned_expected_(2), eps_);
}

// Test NED to ECEF round-trip
TEST_F(AccelerationTransformTest, NedToEcefAcceleration) {
  Eigen::Vector3d a_ecef_result = nt::ned2ecefa(a_ned_expected_, lla0_);
  EXPECT_NEAR(a_ecef_result(0), a_ecef_expected_(0), eps_);
  EXPECT_NEAR(a_ecef_result(1), a_ecef_expected_(1), eps_);
  EXPECT_NEAR(a_ecef_result(2), a_ecef_expected_(2), eps_);
}

// Test ECEF to ENU conversion
TEST_F(AccelerationTransformTest, EcefToEnuAcceleration) {
  Eigen::Vector3d a_enu_result = nt::ecef2enua(a_ecef_expected_, lla0_);
  EXPECT_NEAR(a_enu_result(0), a_enu_expected_(0), eps_);
  EXPECT_NEAR(a_enu_result(1), a_enu_expected_(1), eps_);
  EXPECT_NEAR(a_enu_result(2), a_enu_expected_(2), eps_);
}

// Test ENU to ECEF round-trip
TEST_F(AccelerationTransformTest, EnuToEcefAcceleration) {
  Eigen::Vector3d a_ecef_result = nt::enu2ecefa(a_enu_expected_, lla0_);
  EXPECT_NEAR(a_ecef_result(0), a_ecef_expected_(0), eps_);
  EXPECT_NEAR(a_ecef_result(1), a_ecef_expected_(1), eps_);
  EXPECT_NEAR(a_ecef_result(2), a_ecef_expected_(2), eps_);
}

// Test ECEF to ECI conversion
TEST_F(AccelerationTransformTest, EcefToEciAcceleration) {
  Eigen::Vector3d a_eci_result =
      nt::ecef2ecia(dt_, a_ecef_expected_, v_ecef_expected_, r_ecef_expected_);
  EXPECT_NEAR(a_eci_result(0), a_eci_expected_(0), eps_);
  EXPECT_NEAR(a_eci_result(1), a_eci_expected_(1), eps_);
  EXPECT_NEAR(a_eci_result(2), a_eci_expected_(2), eps_);
}

// Test ECI to ECEF round-trip
TEST_F(AccelerationTransformTest, EciToEcefAcceleration) {
  Eigen::Vector3d r_eci_test = nt::ecef2eci(dt_, r_ecef_expected_);
  Eigen::Vector3d v_eci_test = nt::ecef2eciv(dt_, v_ecef_expected_, r_ecef_expected_);
  Eigen::Vector3d a_ecef_result = nt::eci2ecefa(dt_, a_eci_expected_, v_eci_test, r_eci_test);

  EXPECT_NEAR(a_ecef_result(0), a_ecef_expected_(0), eps_);
  EXPECT_NEAR(a_ecef_result(1), a_ecef_expected_(1), eps_);
  EXPECT_NEAR(a_ecef_result(2), a_ecef_expected_(2), eps_);
}

// Test the full chained conversion from ECI to NED
TEST_F(AccelerationTransformTest, EciToNedAcceleration) {
  Eigen::Vector3d r_eci_test = nt::ecef2eci(dt_, r_ecef_expected_);
  Eigen::Vector3d v_eci_test = nt::ecef2eciv(dt_, v_ecef_expected_, r_ecef_expected_);
  Eigen::Vector3d a_ned_result = nt::eci2neda(dt_, a_eci_expected_, v_eci_test, r_eci_test, lla0_);

  EXPECT_NEAR(a_ned_result(0), a_ned_expected_(0), eps_);
  EXPECT_NEAR(a_ned_result(1), a_ned_expected_(1), eps_);
  EXPECT_NEAR(a_ned_result(2), a_ned_expected_(2), eps_);
}

// Test the full chained conversion from NED to ECI
TEST_F(AccelerationTransformTest, NedToEciAcceleration) {
  Eigen::Vector3d r_ned_test = nt::ecef2ned(r_ecef_expected_, lla0_);
  Eigen::Vector3d v_ned_test = nt::ecef2nedv(v_ecef_expected_, lla0_);
  Eigen::Vector3d a_eci_result = nt::ned2ecia(dt_, a_ned_expected_, v_ned_test, r_ned_test, lla0_);

  EXPECT_NEAR(a_eci_result(0), a_eci_expected_(0), eps_);
  EXPECT_NEAR(a_eci_result(1), a_eci_expected_(1), eps_);
  EXPECT_NEAR(a_eci_result(2), a_eci_expected_(2), eps_);
}

// Test the full chained conversion from ECI to ENU
TEST_F(AccelerationTransformTest, EciToEnuAcceleration) {
  Eigen::Vector3d r_eci_test = nt::ecef2eci(dt_, r_ecef_expected_);
  Eigen::Vector3d v_eci_test = nt::ecef2eciv(dt_, v_ecef_expected_, r_ecef_expected_);
  Eigen::Vector3d a_enu_result = nt::eci2enua(dt_, a_eci_expected_, v_eci_test, r_eci_test, lla0_);

  EXPECT_NEAR(a_enu_result(0), a_enu_expected_(0), eps_);
  EXPECT_NEAR(a_enu_result(1), a_enu_expected_(1), eps_);
  EXPECT_NEAR(a_enu_result(2), a_enu_expected_(2), eps_);
}

// Test the full chained conversion from ENU to ECI
TEST_F(AccelerationTransformTest, EnuToEciAcceleration) {
  Eigen::Vector3d r_enu_test = nt::ecef2enu(r_ecef_expected_, lla0_);
  Eigen::Vector3d v_enu_test = nt::ecef2enuv(v_ecef_expected_, lla0_);
  Eigen::Vector3d a_eci_result = nt::enu2ecia(dt_, a_enu_expected_, v_enu_test, r_enu_test, lla0_);

  EXPECT_NEAR(a_eci_result(0), a_eci_expected_(0), eps_);
  EXPECT_NEAR(a_eci_result(1), a_eci_expected_(1), eps_);
  EXPECT_NEAR(a_eci_result(2), a_eci_expected_(2), eps_);
}
