#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <limits>
#include <navtools/core/attitude-transforms.hpp>
#include <navtools/core/constants.hpp>
#include <random>

// Use a test fixture for shared setup.
class AttitudeTransformsTest : public ::testing::Test {
 protected:
  const double TOL = 10.0 * std::numeric_limits<double>::epsilon();
  std::mt19937_64 rng;
  std::uniform_real_distribution<double> uniform;

  void SetUp() override {
    std::random_device rd;
    rng.seed(rd());
    uniform = std::uniform_real_distribution<double>(-nt::PI<double> + TOL, nt::PI<double> - TOL);
  }

  template <nt::RotationOrder order>
  Eigen::Vector3d GetRandomAngles() {
    if constexpr ((order == nt::RotationOrder::YXZ) || (order == nt::RotationOrder::ZXY)) {
      return Eigen::Vector3d(0.5 * uniform(rng), uniform(rng), uniform(rng));
    } else if constexpr ((order == nt::RotationOrder::XYZ) || (order == nt::RotationOrder::ZYX)) {
      return Eigen::Vector3d(uniform(rng), 0.5 * uniform(rng), uniform(rng));
    } else if constexpr ((order == nt::RotationOrder::XZY) || (order == nt::RotationOrder::YZX)) {
      return Eigen::Vector3d(uniform(rng), uniform(rng), 0.5 * uniform(rng));
    } else if constexpr (
        (order == nt::RotationOrder::XYX) || (order == nt::RotationOrder::XZX) ||
        (order == nt::RotationOrder::YXY) || (order == nt::RotationOrder::YZY) ||
        (order == nt::RotationOrder::ZXZ) || (order == nt::RotationOrder::ZYZ)) {
      return Eigen::Vector3d(uniform(rng), std::abs(0.5 * uniform(rng)), uniform(rng));
    }
    return Eigen::Vector3d::Zero();
  }

  template <nt::RotationOrder order>
  bool CheckGimbalLock(Eigen::Vector3d& truth) {
    if constexpr ((order == nt::RotationOrder::YXZ) || (order == nt::RotationOrder::ZXY)) {
      if (std::abs(std::abs(truth(0)) - nt::HALF_PI<double>) > TOL) {
        return false;
      }
    } else if constexpr (
        (order == nt::RotationOrder::XYZ) || (order == nt::RotationOrder::ZYX) ||
        (order == nt::RotationOrder::XYX) || (order == nt::RotationOrder::XZX) ||
        (order == nt::RotationOrder::YXY) || (order == nt::RotationOrder::YZY) ||
        (order == nt::RotationOrder::ZXZ) || (order == nt::RotationOrder::ZYZ)) {
      if (std::abs(std::abs(truth(1)) - nt::HALF_PI<double>) > TOL) {
        return false;
      }
    } else if constexpr ((order == nt::RotationOrder::XZY) || (order == nt::RotationOrder::YZX)) {
      if (std::abs(std::abs(truth(2)) - nt::HALF_PI<double>) > TOL) {
        return false;
      }
    }
    return true;
  }

  template <nt::RotationOrder order>
  void TestDcmRoundTrip() {
    for (int i = 0; i < 1; ++i) {
      Eigen::Vector3d true_angles = GetRandomAngles<order>();
      Eigen::Matrix3d R = nt::euler2dcm<order, Eigen::Matrix3d>(true_angles);
      Eigen::Vector3d est_angles = nt::dcm2euler<order, Eigen::Vector3d>(R);

      // Check for gimbal lock, where the sum or difference of angles is preserved.
      if (CheckGimbalLock<order>(true_angles)) {
        EXPECT_NEAR(true_angles(0), est_angles(0), TOL);
        EXPECT_NEAR(true_angles(1), est_angles(1), TOL);
        EXPECT_NEAR(true_angles(2), est_angles(2), TOL);
      }
    }
  }

  template <nt::RotationOrder order>
  void TestQuatRoundTrip() {
    for (int i = 0; i < 100; ++i) {
      Eigen::Vector3d true_angles = GetRandomAngles<order>();
      Eigen::Vector4d Q = nt::euler2quat<order, Eigen::Vector4d>(true_angles);
      Eigen::Vector3d est_angles = nt::quat2euler<order, Eigen::Vector3d>(Q);

      // Check for gimbal lock, where the sum or difference of angles is preserved.
      if (CheckGimbalLock<order>(true_angles)) {
        EXPECT_NEAR(true_angles(0), est_angles(0), TOL);
        EXPECT_NEAR(true_angles(1), est_angles(1), TOL);
        EXPECT_NEAR(true_angles(2), est_angles(2), TOL);
      }
    }
  }

  template <nt::RotationOrder order>
  void TestQuat2DcmTrip() {
    for (int i = 0; i < 100; ++i) {
      Eigen::Vector3d true_angles = GetRandomAngles<order>();
      Eigen::Vector4d Q = nt::euler2quat<order, Eigen::Vector4d>(true_angles);
      Eigen::Matrix3d R = nt::quat2dcm<Eigen::Matrix3d>(Q);
      Eigen::Vector3d est_angles = nt::dcm2euler<order, Eigen::Vector3d>(R);

      // Check for gimbal lock, where the sum or difference of angles is preserved.
      if (CheckGimbalLock<order>(true_angles)) {
        EXPECT_NEAR(true_angles(0), est_angles(0), TOL);
        EXPECT_NEAR(true_angles(1), est_angles(1), TOL);
        EXPECT_NEAR(true_angles(2), est_angles(2), TOL);
      }
    }
  }

  template <nt::RotationOrder order>
  void TestDcm2QuatTrip() {
    for (int i = 0; i < 100; ++i) {
      Eigen::Vector3d true_angles = GetRandomAngles<order>();
      Eigen::Matrix3d R = nt::euler2dcm<order, Eigen::Matrix3d>(true_angles);
      Eigen::Vector4d Q = nt::dcm2quat<Eigen::Vector4d>(R);
      Eigen::Vector3d est_angles = nt::quat2euler<order, Eigen::Vector3d>(Q);

      // Check for gimbal lock, where the sum or difference of angles is preserved.
      if (CheckGimbalLock<order>(true_angles)) {
        EXPECT_NEAR(true_angles(0), est_angles(0), TOL);
        EXPECT_NEAR(true_angles(1), est_angles(1), TOL);
        EXPECT_NEAR(true_angles(2), est_angles(2), TOL);
      }
    }
  }
};

// Test all 12 rotation orders for DCM conversions.
TEST_F(AttitudeTransformsTest, All_Dcm_Rotation_Orders_Roundtrip) {
  TestDcmRoundTrip<nt::RotationOrder::XYZ>();
  TestDcmRoundTrip<nt::RotationOrder::XZY>();
  TestDcmRoundTrip<nt::RotationOrder::YXZ>();
  TestDcmRoundTrip<nt::RotationOrder::YZX>();
  TestDcmRoundTrip<nt::RotationOrder::ZXY>();
  TestDcmRoundTrip<nt::RotationOrder::ZYX>();
  TestDcmRoundTrip<nt::RotationOrder::XYX>();
  TestDcmRoundTrip<nt::RotationOrder::XZX>();
  TestDcmRoundTrip<nt::RotationOrder::YXY>();
  TestDcmRoundTrip<nt::RotationOrder::YZY>();
  TestDcmRoundTrip<nt::RotationOrder::ZXZ>();
  TestDcmRoundTrip<nt::RotationOrder::ZYZ>();
}

// Test all 12 rotation orders for Quaternion conversions.
TEST_F(AttitudeTransformsTest, All_Quat_Rotation_Orders_Roundtrip) {
  TestQuatRoundTrip<nt::RotationOrder::ZYX>();
  TestQuatRoundTrip<nt::RotationOrder::ZXY>();
  TestQuatRoundTrip<nt::RotationOrder::YZX>();
  TestQuatRoundTrip<nt::RotationOrder::YXZ>();
  TestQuatRoundTrip<nt::RotationOrder::XYZ>();
  TestQuatRoundTrip<nt::RotationOrder::XZY>();
  TestQuatRoundTrip<nt::RotationOrder::ZYZ>();
  TestQuatRoundTrip<nt::RotationOrder::ZXZ>();
  TestQuatRoundTrip<nt::RotationOrder::YZY>();
  TestQuatRoundTrip<nt::RotationOrder::YXY>();
  TestQuatRoundTrip<nt::RotationOrder::XZX>();
  TestQuatRoundTrip<nt::RotationOrder::XYX>();
}

// Test dcm2quat
TEST_F(AttitudeTransformsTest, All_Dcm_to_Quat_Transforms) {
  TestDcm2QuatTrip<nt::RotationOrder::XYZ>();
  TestDcm2QuatTrip<nt::RotationOrder::XZY>();
  TestDcm2QuatTrip<nt::RotationOrder::YXZ>();
  TestDcm2QuatTrip<nt::RotationOrder::YZX>();
  TestDcm2QuatTrip<nt::RotationOrder::ZXY>();
  TestDcm2QuatTrip<nt::RotationOrder::ZYX>();
  TestDcm2QuatTrip<nt::RotationOrder::XYX>();
  TestDcm2QuatTrip<nt::RotationOrder::XZX>();
  TestDcm2QuatTrip<nt::RotationOrder::YXY>();
  TestDcm2QuatTrip<nt::RotationOrder::YZY>();
  TestDcm2QuatTrip<nt::RotationOrder::ZXZ>();
  TestDcm2QuatTrip<nt::RotationOrder::ZYZ>();
}

// Test quat2dcm
TEST_F(AttitudeTransformsTest, All_Quat_to_Dcm_Transforms) {
  TestQuat2DcmTrip<nt::RotationOrder::ZYX>();
  TestQuat2DcmTrip<nt::RotationOrder::ZXY>();
  TestQuat2DcmTrip<nt::RotationOrder::YZX>();
  TestQuat2DcmTrip<nt::RotationOrder::YXZ>();
  TestQuat2DcmTrip<nt::RotationOrder::XYZ>();
  TestQuat2DcmTrip<nt::RotationOrder::XZY>();
  TestQuat2DcmTrip<nt::RotationOrder::ZYZ>();
  TestQuat2DcmTrip<nt::RotationOrder::ZXZ>();
  TestQuat2DcmTrip<nt::RotationOrder::YZY>();
  TestQuat2DcmTrip<nt::RotationOrder::YXY>();
  TestQuat2DcmTrip<nt::RotationOrder::XZX>();
  TestQuat2DcmTrip<nt::RotationOrder::XYX>();
}