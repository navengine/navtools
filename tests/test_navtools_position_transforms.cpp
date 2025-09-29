#include <Eigen/Dense>
#include <cmath>
#include <gtest/gtest.h>
#include <navtools/core/coordinate-dcms.hpp>
#include <navtools/core/position-transforms.hpp>

// Define a test fixture to hold common data
class CoordinateTransformTest : public ::testing::Test {
 protected:
  // Tolerance for floating-point comparisons
  const double kTolerance = 1e-9;

  // Test points in LLA format (lat, lon, alt)
  // Los Angeles, CA (34.0517,	-118.2452,	70.745)
  const Eigen::Vector3d lla_la_ = {0.594314087066431, -2.06376826264174, 70.7450907109305};
  // Oslo, Norway (57.4638, 10.3407, 100.391)
  const Eigen::Vector3d lla_oslo_ = {1.00293278714122, 0.18047974773271, 100.391102797352};
  // North Pole
  const Eigen::Vector3d lla_pole_north_ = {nt::HALF_PI<double>, 0.0, 0.0};
  // On Equator, 1000km up (rad, rad, m)
  const Eigen::Vector3d lla_high_alt_ = {0.0, 0.0, 1000000.0};

  // Test points in ECEF format (x, y, z)
  const Eigen::Vector3d ecef_la_ = {-2503523.000000021, -4660216.9999999879, 3551238.0000000047};
  const Eigen::Vector3d ecef_oslo_ = {3382760.0000000019, 617235.99999999872, 5353941.0000000019};
  const Eigen::Vector3d ecef_origin_ = {0.0, 0.0, 0.0};

  // Test point in NED format
  const Eigen::Vector3d ned_test_ = {100.0, -50.0, 20.0};

  // Test point in ENU format
  const Eigen::Vector3d enu_test_ = {100.0, -50.0, 20.0};
};

// LLA <-> ECEF tests
TEST_F(CoordinateTransformTest, lla2ecefAndBack) {
  Eigen::Vector3d xyz, lla;

  // Test 1: Los Angeles, CA
  nt::lla2ecef(lla_la_, xyz);
  EXPECT_NEAR(xyz(0), ecef_la_(0), kTolerance);
  EXPECT_NEAR(xyz(1), ecef_la_(1), kTolerance);
  EXPECT_NEAR(xyz(2), ecef_la_(2), kTolerance);

  // Round trip
  nt::ecef2lla(xyz, lla);
  EXPECT_NEAR(lla(0), lla_la_(0), kTolerance);
  EXPECT_NEAR(lla(1), lla_la_(1), kTolerance);
  EXPECT_NEAR(lla(2), lla_la_(2), kTolerance);

  // Test 2: Oslo, Norway (High Latitude)
  nt::lla2ecef(lla_oslo_, xyz);
  EXPECT_NEAR(xyz(0), ecef_oslo_(0), kTolerance);
  EXPECT_NEAR(xyz(1), ecef_oslo_(1), kTolerance);
  EXPECT_NEAR(xyz(2), ecef_oslo_(2), kTolerance);

  // Round trip
  nt::ecef2lla(xyz, lla);
  EXPECT_NEAR(lla(0), lla_oslo_(0), kTolerance);
  EXPECT_NEAR(lla(1), lla_oslo_(1), kTolerance);
  EXPECT_NEAR(lla(2), lla_oslo_(2), kTolerance);

  // Test 3: North Pole
  nt::lla2ecef(lla_pole_north_, xyz);
  nt::ecef2lla(xyz, lla);
  EXPECT_NEAR(lla(0), lla_pole_north_(0), kTolerance);
  // Longitude is undefined at the poles, so skip the check
  EXPECT_NEAR(lla(2), lla_pole_north_(2), kTolerance);
}

// ECEF <-> ECI tests
TEST_F(CoordinateTransformTest, ecef2eciAndBack) {
  // Test time elapsed (e.g., 6 hours)
  const double dt = 6 * 3600.0;
  Eigen::Vector3d eci_result, ecef_round_trip;

  // ECEF to ECI
  nt::ecef2eci(dt, ecef_la_, eci_result);

  // ECI to ECEF (round trip)
  nt::eci2ecef(dt, eci_result, ecef_round_trip);
  EXPECT_NEAR(ecef_round_trip(0), ecef_la_(0), kTolerance);
  EXPECT_NEAR(ecef_round_trip(1), ecef_la_(1), kTolerance);
  EXPECT_NEAR(ecef_round_trip(2), ecef_la_(2), kTolerance);
}

// Checks the rotation between ECEF frames over time via the ECI frame
TEST_F(CoordinateTransformTest, ecefRotationViaEciCheck) {
  // Use a convenient time interval, e.g., 6 hours
  const double dt1 = 6.0 * 3600.0;
  const double dt2 = 12.0 * 3600.0;

  Eigen::Vector3d eci_result, ecef_result;

  // 1. Convert initial ECEF point to ECI at time dt1
  nt::ecef2eci(dt1, ecef_la_, eci_result);

  // 2. Convert the ECI point back to ECEF at a later time, dt2
  nt::eci2ecef(dt2, eci_result, ecef_result);

  // 3. Manually calculate the expected final ECEF vector
  // The final ECEF vector should be the initial one rotated by the time difference (dt2 - dt1)
  const double rotation_angle = -nt::WGS84_OMEGA<double> * (dt2 - dt1);
  Eigen::Matrix3d expected_dcm;
  expected_dcm << std::cos(rotation_angle), -std::sin(rotation_angle), 0, std::sin(rotation_angle),
      std::cos(rotation_angle), 0, 0, 0, 1;
  Eigen::Vector3d expected_ecef = expected_dcm * ecef_la_;

  // 4. Assert that the calculated ECEF matches the expected one
  EXPECT_NEAR(ecef_result(0), expected_ecef(0), kTolerance);
  EXPECT_NEAR(ecef_result(1), expected_ecef(1), kTolerance);
  EXPECT_NEAR(ecef_result(2), expected_ecef(2), kTolerance);
}

// LLA <-> ECI tests
TEST_F(CoordinateTransformTest, lla2eciAndBack) {
  const double dt = 12 * 3600.0;
  Eigen::Vector3d eci_result, lla_round_trip;

  // LLA to ECI
  nt::lla2eci(dt, lla_la_, eci_result);

  // ECI to LLA (round trip)
  nt::eci2lla(dt, eci_result, lla_round_trip);

  // Check that the final LLA is the same as the initial one
  EXPECT_NEAR(lla_round_trip(0), lla_la_(0), kTolerance);
  EXPECT_NEAR(lla_round_trip(1), lla_la_(1), kTolerance);
  EXPECT_NEAR(lla_round_trip(2), lla_la_(2), kTolerance);
}

// ECEF <-> NED tests
TEST_F(CoordinateTransformTest, ecef2nedAndBack) {
  // Use a known ECEF point and a reference LLA point
  Eigen::Vector3d lla0 = lla_la_;
  Eigen::Vector3d ned_result, ecef_round_trip;

  // ECEF -> LLA0
  Eigen::Vector3d ecef0;
  nt::lla2ecef(lla0, ecef0);

  // Compute a new ECEF point offset by the NED vector
  Eigen::Vector3d new_ecef;
  nt::ned2ecef(ned_test_, lla0, new_ecef);

  // ECEF -> NED (round trip)
  nt::ecef2ned(new_ecef, lla0, ned_result);
  EXPECT_NEAR(ned_result(0), ned_test_(0), kTolerance);
  EXPECT_NEAR(ned_result(1), ned_test_(1), kTolerance);
  EXPECT_NEAR(ned_result(2), ned_test_(2), kTolerance);
}

// ECEF <-> ENU tests
TEST_F(CoordinateTransformTest, ecef2enuAndBack) {
  // Use a known ECEF point and a reference LLA point
  Eigen::Vector3d lla0 = lla_la_;
  Eigen::Vector3d enu_result, ecef_round_trip;

  // Compute a new ECEF point offset by the ENU vector
  Eigen::Vector3d new_ecef;
  nt::enu2ecef(enu_test_, lla0, new_ecef);

  // ECEF -> ENU (round trip)
  nt::ecef2enu(new_ecef, lla0, enu_result);
  EXPECT_NEAR(enu_result(0), enu_test_(0), kTolerance);
  EXPECT_NEAR(enu_result(1), enu_test_(1), kTolerance);
  EXPECT_NEAR(enu_result(2), enu_test_(2), kTolerance);
}

// NED <-> ENU tests
TEST_F(CoordinateTransformTest, ned2enuAndBack) {
  Eigen::Vector3d enu_result, ned_round_trip;

  // NED -> ENU
  nt::ned2enu(ned_test_, enu_result);

  // ENU -> NED (round trip)
  nt::enu2ned(enu_result, ned_round_trip);
  EXPECT_NEAR(ned_round_trip(0), ned_test_(0), kTolerance);
  EXPECT_NEAR(ned_round_trip(1), ned_test_(1), kTolerance);
  EXPECT_NEAR(ned_round_trip(2), ned_test_(2), kTolerance);
}

// Comprehensive chain tests
TEST_F(CoordinateTransformTest, lla2nedViaEcefAndBack) {
  // Chain: lla -> ecef -> ned -> ecef -> lla
  Eigen::Vector3d lla0 = lla_la_;
  Eigen::Vector3d xyz, ned, xyz_roundtrip, lla_roundtrip;

  // lla to ecef
  nt::lla2ecef(lla0, xyz);

  // offset by a local ned vector
  nt::ned2ecef(ned_test_, lla0, xyz_roundtrip);

  // convert the new ecef to ned
  nt::ecef2ned(xyz_roundtrip, lla0, ned);

  // Verify the NED result
  EXPECT_NEAR(ned(0), ned_test_(0), kTolerance);
  EXPECT_NEAR(ned(1), ned_test_(1), kTolerance);
  EXPECT_NEAR(ned(2), ned_test_(2), kTolerance);

  // And finally back to LLA to ensure the conversion is correct
  nt::ecef2lla(xyz_roundtrip, lla_roundtrip);
}