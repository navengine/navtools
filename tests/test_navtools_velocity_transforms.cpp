#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <navtools/core/coordinate-dcms.hpp>
#include <navtools/core/position-transforms.hpp>
#include <navtools/core/velocity-transforms.hpp>

// Define a test fixture to hold common test data
class VelocityTransformTest : public ::testing::Test {
 protected:
  // Tolerance for floating-point comparisons
  const double eps_ = 1e-9;

  // Time elapsed for ECI transforms
  const double dt_ = 3600.0;  // 1 hour

  // Reference point (Los Angeles, CA)
  const Eigen::Vector3d lla_la_ = {0.594314087066431, -2.06376826264174, 70.745};
  const Eigen::Vector3d lla0_ = {0.594, -2.063, 0.0};

  // Reference position vectors
  Eigen::Vector3d r_ecef_expected_ = nt::lla2ecef(lla_la_);
  Eigen::Vector3d r_eci_expected_ = nt::lla2eci(dt_, lla_la_);
  Eigen::Vector3d r_ned_expected_ = nt::lla2ned(lla_la_, lla0_);
  Eigen::Vector3d r_enu_expected_ = nt::lla2enu(lla_la_, lla0_);

  // A velocity vector to test with (m/s)
  // Non-zero components to fully test all calculations
  Eigen::Vector3d v_ecef_expected_ = {1.23, -4.56, 7.89};

  // Pre-calculated values for comparison (based on the formulas)
  Eigen::Vector3d v_ned_expected_ = nt::ecef2nedDcm<Eigen::Matrix3d>(lla0_) * v_ecef_expected_;
  Eigen::Vector3d v_enu_expected_ = nt::ecef2enuDcm<Eigen::Matrix3d>(lla0_) * v_ecef_expected_;
  Eigen::Vector3d v_eci_expected_ =
      nt::ecef2eciDcm<Eigen::Matrix3d>(dt_) *
      (v_ecef_expected_ + nt::OMEGA_ECEF<double>.cross(r_ecef_expected_));
};

// Test ECEF to NED conversion
TEST_F(VelocityTransformTest, ecefToNedVelocity) {
  Eigen::Vector3d v_ned_result = nt::ecef2nedv(v_ecef_expected_, lla0_);
  EXPECT_NEAR(v_ned_result(0), v_ned_expected_(0), eps_);
  EXPECT_NEAR(v_ned_result(1), v_ned_expected_(1), eps_);
  EXPECT_NEAR(v_ned_result(2), v_ned_expected_(2), eps_);
}

// Test NED to ECEF round-trip
TEST_F(VelocityTransformTest, nedToEcefVelocity) {
  Eigen::Vector3d v_ecef_result = nt::ned2ecefv(v_ned_expected_, lla0_);
  EXPECT_NEAR(v_ecef_result(0), v_ecef_expected_(0), eps_);
  EXPECT_NEAR(v_ecef_result(1), v_ecef_expected_(1), eps_);
  EXPECT_NEAR(v_ecef_result(2), v_ecef_expected_(2), eps_);
}

// Test ECEF to ENU conversion
TEST_F(VelocityTransformTest, ecefToEnuVelocity) {
  Eigen::Vector3d v_enu_result = nt::ecef2enuv(v_ecef_expected_, lla0_);
  EXPECT_NEAR(v_enu_result(0), v_enu_expected_(0), eps_);
  EXPECT_NEAR(v_enu_result(1), v_enu_expected_(1), eps_);
  EXPECT_NEAR(v_enu_result(2), v_enu_expected_(2), eps_);
}

// Test ENU to ECEF round-trip
TEST_F(VelocityTransformTest, enuToEcefVelocity) {
  Eigen::Vector3d v_ecef_result = nt::enu2ecefv(v_enu_expected_, lla0_);
  EXPECT_NEAR(v_ecef_result(0), v_ecef_expected_(0), eps_);
  EXPECT_NEAR(v_ecef_result(1), v_ecef_expected_(1), eps_);
  EXPECT_NEAR(v_ecef_result(2), v_ecef_expected_(2), eps_);
}

// Test ECEF to ECI conversion
TEST_F(VelocityTransformTest, ecefToEciVelocity) {
  Eigen::Vector3d v_eci_result = nt::ecef2eciv(dt_, v_ecef_expected_, r_ecef_expected_);
  EXPECT_NEAR(v_eci_result(0), v_eci_expected_(0), eps_);
  EXPECT_NEAR(v_eci_result(1), v_eci_expected_(1), eps_);
  EXPECT_NEAR(v_eci_result(2), v_eci_expected_(2), eps_);
}

// Test ECI to ECEF round-trip
TEST_F(VelocityTransformTest, eciToEcefVelocity) {
  Eigen::Vector3d v_ecef_result = nt::eci2ecefv(dt_, v_eci_expected_, r_eci_expected_);
  EXPECT_NEAR(v_ecef_result(0), v_ecef_expected_(0), eps_);
  EXPECT_NEAR(v_ecef_result(1), v_ecef_expected_(1), eps_);
  EXPECT_NEAR(v_ecef_result(2), v_ecef_expected_(2), eps_);
}

// Test the full chained conversion from ECI to NED
TEST_F(VelocityTransformTest, eciToNedVelocity) {
  Eigen::Vector3d v_ned_result = nt::eci2nedv(dt_, v_eci_expected_, r_eci_expected_, lla0_);
  EXPECT_NEAR(v_ned_result(0), v_ned_expected_(0), eps_);
  EXPECT_NEAR(v_ned_result(1), v_ned_expected_(1), eps_);
  EXPECT_NEAR(v_ned_result(2), v_ned_expected_(2), eps_);
}

// Test the full chained conversion from NED to ECI
TEST_F(VelocityTransformTest, nedToEciVelocity) {
  Eigen::Vector3d v_eci_result = nt::ned2eciv(dt_, v_ned_expected_, r_ned_expected_, lla0_);
  EXPECT_NEAR(v_eci_result(0), v_eci_expected_(0), eps_);
  EXPECT_NEAR(v_eci_result(1), v_eci_expected_(1), eps_);
  EXPECT_NEAR(v_eci_result(2), v_eci_expected_(2), eps_);
}

// Test the full chained conversion from ECI to ENU
TEST_F(VelocityTransformTest, eciToEnuVelocity) {
  Eigen::Vector3d v_enu_result = nt::eci2enuv(dt_, v_eci_expected_, r_eci_expected_, lla0_);
  EXPECT_NEAR(v_enu_result(0), v_enu_expected_(0), eps_);
  EXPECT_NEAR(v_enu_result(1), v_enu_expected_(1), eps_);
  EXPECT_NEAR(v_enu_result(2), v_enu_expected_(2), eps_);
}

// Test the full chained conversion from NED to ENU
TEST_F(VelocityTransformTest, enuToEciVelocity) {
  Eigen::Vector3d v_eci_result = nt::enu2eciv(dt_, v_enu_expected_, r_enu_expected_, lla0_);
  EXPECT_NEAR(v_eci_result(0), v_eci_expected_(0), eps_);
  EXPECT_NEAR(v_eci_result(1), v_eci_expected_(1), eps_);
  EXPECT_NEAR(v_eci_result(2), v_eci_expected_(2), eps_);
}
