#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <navtools/core/acceleration-transforms.hpp>
#include <navtools/core/angular-rate-transforms.hpp>
#include <navtools/core/coordinate-dcms.hpp>

// A test fixture to hold common data for all tests.
class AngularVelocityTest : public ::testing::Test {
 protected:
  // Tolerance for floating-point comparisons
  const double eps_ = 1e-9;

  // Time elapsed for ECI transforms
  const double dt_ = 3600.0;  // 1 hour

  // Reference point (Los Angeles, CA)
  const Eigen::Vector3d lla_la_ = {0.594314087066431, -2.06376826264174, 70.745};
  const Eigen::Vector3d lla0_ = {0.594, -2.063, 0.0};

  Eigen::Vector3d w_ecef_expected_ = {1.0, -0.2, 0.03};
  Eigen::Vector3d w_ned_expected_ = nt::ecef2nedDcm(lla0_) * w_ecef_expected_;
  Eigen::Vector3d w_enu_expected_ = nt::ecef2enuDcm(lla0_) * w_ecef_expected_;
  Eigen::Vector3d w_eci_expected_ =
      nt::ecef2eciDcm(dt_) * (w_ecef_expected_ + nt::OMEGA_ECEF<double>);
};

// Test for ecef2eciw
TEST_F(AngularVelocityTest, Ecef2EciAngularRate) {
  Eigen::Vector3d w_eci_actual = nt::ecef2eciw(dt_, w_ecef_expected_);
  EXPECT_NEAR(w_eci_actual(0), w_eci_expected_(0), eps_);
  EXPECT_NEAR(w_eci_actual(1), w_eci_expected_(1), eps_);
  EXPECT_NEAR(w_eci_actual(2), w_eci_expected_(2), eps_);
}

// Test for eci2ecefw
TEST_F(AngularVelocityTest, Eci2EcefAngularRate) {
  Eigen::Vector3d w_xyz_actual = nt::eci2ecefw(dt_, w_eci_expected_);
  EXPECT_NEAR(w_xyz_actual(0), w_ecef_expected_(0), eps_);
  EXPECT_NEAR(w_xyz_actual(1), w_ecef_expected_(1), eps_);
  EXPECT_NEAR(w_xyz_actual(2), w_ecef_expected_(2), eps_);
}

// Test for ecef2nedw
TEST_F(AngularVelocityTest, Ecef2NedAngularRate) {
  Eigen::Vector3d w_ned_actual = nt::ecef2nedw(w_ecef_expected_, lla0_);
  EXPECT_NEAR(w_ned_actual(0), w_ned_expected_(0), eps_);
  EXPECT_NEAR(w_ned_actual(1), w_ned_expected_(1), eps_);
  EXPECT_NEAR(w_ned_actual(2), w_ned_expected_(2), eps_);
}

// Test for ecef2enuw
TEST_F(AngularVelocityTest, Ecef2EnuAngularRate) {
  Eigen::Vector3d w_enu_actual = nt::ecef2enuw(w_ecef_expected_, lla0_);
  EXPECT_NEAR(w_enu_actual(0), w_enu_expected_(0), eps_);
  EXPECT_NEAR(w_enu_actual(1), w_enu_expected_(1), eps_);
  EXPECT_NEAR(w_enu_actual(2), w_enu_expected_(2), eps_);
}

// Test for ned2ecefw
TEST_F(AngularVelocityTest, Ned2EcefAngularRate) {
  Eigen::Vector3d w_xyz_actual = nt::ned2ecefw(w_ned_expected_, lla0_);
  EXPECT_NEAR(w_xyz_actual(0), w_ecef_expected_(0), eps_);
  EXPECT_NEAR(w_xyz_actual(1), w_ecef_expected_(1), eps_);
  EXPECT_NEAR(w_xyz_actual(2), w_ecef_expected_(2), eps_);
}

// Test for enu2ecefw
TEST_F(AngularVelocityTest, Enu2EcefAngularRate) {
  Eigen::Vector3d w_xyz_actual = nt::enu2ecefw(w_enu_expected_, lla0_);
  EXPECT_NEAR(w_xyz_actual(0), w_ecef_expected_(0), eps_);
  EXPECT_NEAR(w_xyz_actual(1), w_ecef_expected_(1), eps_);
  EXPECT_NEAR(w_xyz_actual(2), w_ecef_expected_(2), eps_);
}

// Test for eci2nedw
TEST_F(AngularVelocityTest, Eci2NedAngularRate) {
  Eigen::Vector3d w_ned_actual = nt::eci2nedw(dt_, w_eci_expected_, lla0_);
  EXPECT_NEAR(w_ned_actual(0), w_ned_expected_(0), eps_);
  EXPECT_NEAR(w_ned_actual(1), w_ned_expected_(1), eps_);
  EXPECT_NEAR(w_ned_actual(2), w_ned_expected_(2), eps_);
}

// Test for eci2enuw
TEST_F(AngularVelocityTest, Eci2EnuAngularRate) {
  Eigen::Vector3d w_enu_actual = nt::eci2enuw(dt_, w_eci_expected_, lla0_);
  EXPECT_NEAR(w_enu_actual(0), w_enu_expected_(0), eps_);
  EXPECT_NEAR(w_enu_actual(1), w_enu_expected_(1), eps_);
  EXPECT_NEAR(w_enu_actual(2), w_enu_expected_(2), eps_);
}

// Test for ned2eciw
TEST_F(AngularVelocityTest, Ned2EciAngularRate) {
  Eigen::Vector3d w_eci_actual = nt::ned2eciw(dt_, w_ned_expected_, lla0_);
  EXPECT_NEAR(w_eci_actual(0), w_eci_expected_(0), eps_);
  EXPECT_NEAR(w_eci_actual(1), w_eci_expected_(1), eps_);
  EXPECT_NEAR(w_eci_actual(2), w_eci_expected_(2), eps_);
}

// Test for enu2eciw
TEST_F(AngularVelocityTest, Enu2EciAngularRate) {
  Eigen::Vector3d w_eci_actual = nt::enu2eciw(dt_, w_enu_expected_, lla0_);
  EXPECT_NEAR(w_eci_actual(0), w_eci_expected_(0), eps_);
  EXPECT_NEAR(w_eci_actual(1), w_eci_expected_(1), eps_);
  EXPECT_NEAR(w_eci_actual(2), w_eci_expected_(2), eps_);
}
