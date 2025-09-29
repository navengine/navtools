#include <gtest/gtest.h>
#include <navtools/attitude>

const double TOL = 10.0 * std::numeric_limits<double>::epsilon();

TEST(ConstructorTest, EulerAnglesTest) {
  nt::EulerAngles<double, nt::RotationOrder::YXZ> eul_empty;
  EXPECT_NEAR(eul_empty(0), 0.0, TOL);
  EXPECT_NEAR(eul_empty(1), 0.0, TOL);
  EXPECT_NEAR(eul_empty(2), 0.0, TOL);

  nt::EulerAngles<double, nt::RotationOrder::ZXY> eul_from_angle(1.0, 2.0, 3.0);
  EXPECT_NEAR(eul_from_angle(0), 1.0, TOL);
  EXPECT_NEAR(eul_from_angle(1), 2.0, TOL);
  EXPECT_NEAR(eul_from_angle(2), 3.0, TOL);

  Eigen::Vector3d vec(0.1, -0.2, 0.3);
  nt::EulerAngles<double, nt::RotationOrder::XZY> eul_from_vec(vec);
  EXPECT_NEAR(eul_from_vec(0), 0.1, TOL);
  EXPECT_NEAR(eul_from_vec(1), -0.2, TOL);
  EXPECT_NEAR(eul_from_vec(2), 0.3, TOL);

  nt::RotationMatrix<double> dcm;
  nt::angles2dcm<nt::RotationOrder::XYZ>(
      -nt::PI<double> / 7.0, nt::PI<double> / 5.0, -nt::PI<double> / 2.0, dcm);
  nt::EulerAngles<double, nt::RotationOrder::XYZ> eul_from_dcm(dcm);
  EXPECT_NEAR(eul_from_dcm(0), -nt::PI<double> / 7.0, TOL);
  EXPECT_NEAR(eul_from_dcm(1), nt::PI<double> / 5.0, TOL);
  EXPECT_NEAR(eul_from_dcm(2), -nt::PI<double> / 2.0, TOL);

  nt::Quaternion<double> quat;
  nt::angles2quat<nt::RotationOrder::ZYX>(
      nt::PI<double> / 4.0, -nt::PI<double> / 6.0, nt::PI<double> / 3.0, quat);
  nt::EulerAngles<double, nt::RotationOrder::ZYX> eul_from_quat(quat);
  EXPECT_NEAR(eul_from_quat(0), nt::PI<double> / 4.0, TOL);
  EXPECT_NEAR(eul_from_quat(1), -nt::PI<double> / 6.0, TOL);
  EXPECT_NEAR(eul_from_quat(2), nt::PI<double> / 3.0, TOL);
}

TEST(ConstructorTest, RotationMatrixTest) {
  nt::RotationMatrix<double> dcm_empty;
  EXPECT_NEAR(dcm_empty(0, 0), 0.0, TOL);
  EXPECT_NEAR(dcm_empty(0, 1), 0.0, TOL);
  EXPECT_NEAR(dcm_empty(0, 2), 0.0, TOL);
  EXPECT_NEAR(dcm_empty(1, 0), 0.0, TOL);
  EXPECT_NEAR(dcm_empty(1, 1), 0.0, TOL);
  EXPECT_NEAR(dcm_empty(1, 2), 0.0, TOL);
  EXPECT_NEAR(dcm_empty(2, 0), 0.0, TOL);
  EXPECT_NEAR(dcm_empty(2, 1), 0.0, TOL);
  EXPECT_NEAR(dcm_empty(2, 2), 0.0, TOL);

  Eigen::Matrix3d mat;
  mat << 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
  nt::RotationMatrix<double> dcm_from_mat(mat);
  EXPECT_NEAR(dcm_from_mat(0, 0), 0.0, TOL);
  EXPECT_NEAR(dcm_from_mat(0, 1), 1.0, TOL);
  EXPECT_NEAR(dcm_from_mat(0, 2), 0.0, TOL);
  EXPECT_NEAR(dcm_from_mat(1, 0), 1.0, TOL);
  EXPECT_NEAR(dcm_from_mat(1, 1), 0.0, TOL);
  EXPECT_NEAR(dcm_from_mat(1, 2), 0.0, TOL);
  EXPECT_NEAR(dcm_from_mat(2, 0), 0.0, TOL);
  EXPECT_NEAR(dcm_from_mat(2, 1), 0.0, TOL);
  EXPECT_NEAR(dcm_from_mat(2, 2), -1.0, TOL);

  nt::EulerAngles<double, nt::RotationOrder::ZYX> eul(
      -nt::PI<double> / 7.0, nt::PI<double> / 5.0, -nt::PI<double> / 2.0);
  Eigen::Matrix3d true_dcm;
  nt::euler2dcm<nt::RotationOrder::ZYX>(eul, true_dcm);
  nt::RotationMatrix<double> dcm_from_eul(eul);
  EXPECT_NEAR(dcm_from_eul(0, 0), true_dcm(0, 0), TOL);
  EXPECT_NEAR(dcm_from_eul(0, 1), true_dcm(0, 1), TOL);
  EXPECT_NEAR(dcm_from_eul(0, 2), true_dcm(0, 2), TOL);
  EXPECT_NEAR(dcm_from_eul(1, 0), true_dcm(1, 0), TOL);
  EXPECT_NEAR(dcm_from_eul(1, 1), true_dcm(1, 1), TOL);
  EXPECT_NEAR(dcm_from_eul(1, 2), true_dcm(1, 2), TOL);
  EXPECT_NEAR(dcm_from_eul(2, 0), true_dcm(2, 0), TOL);
  EXPECT_NEAR(dcm_from_eul(2, 1), true_dcm(2, 1), TOL);
  EXPECT_NEAR(dcm_from_eul(2, 2), true_dcm(2, 2), TOL);

  nt::Quaternion<double> quat(0.0, 0.0, 0.0, 1.0);
  nt::RotationMatrix<double> dcm_from_quat(quat);
  EXPECT_NEAR(dcm_from_quat(0, 0), 1.0, TOL);
  EXPECT_NEAR(dcm_from_quat(0, 1), 0.0, TOL);
  EXPECT_NEAR(dcm_from_quat(0, 2), 0.0, TOL);
  EXPECT_NEAR(dcm_from_quat(1, 0), 0.0, TOL);
  EXPECT_NEAR(dcm_from_quat(1, 1), 1.0, TOL);
  EXPECT_NEAR(dcm_from_quat(1, 2), 0.0, TOL);
  EXPECT_NEAR(dcm_from_quat(2, 0), 0.0, TOL);
  EXPECT_NEAR(dcm_from_quat(2, 1), 0.0, TOL);
  EXPECT_NEAR(dcm_from_quat(2, 2), 1.0, TOL);
}

TEST(ConstructorTest, QuaternionTest) {
  nt::Quaternion<double> quat_empty;
  EXPECT_NEAR(quat_empty(0), 0.0, TOL);
  EXPECT_NEAR(quat_empty(1), 0.0, TOL);
  EXPECT_NEAR(quat_empty(2), 0.0, TOL);
  EXPECT_NEAR(quat_empty(3), 0.0, TOL);

  Eigen::Vector4d vec(-0.1, 0.2, -0.3, 0.4);
  vec /= vec.norm();
  nt::Quaternion<double> quat_from_vec(vec);
  EXPECT_NEAR(quat_from_vec(0), vec(0), TOL);
  EXPECT_NEAR(quat_from_vec(1), vec(1), TOL);
  EXPECT_NEAR(quat_from_vec(2), vec(2), TOL);
  EXPECT_NEAR(quat_from_vec(3), vec(3), TOL);

  nt::EulerAngles<double, nt::RotationOrder::ZYX> eul(
      -nt::PI<double> / 7.0, nt::PI<double> / 5.0, -nt::PI<double> / 2.0);
  Eigen::Vector4d true_quat;
  nt::euler2quat<nt::RotationOrder::ZYX>(eul, true_quat);
  nt::Quaternion<double> quat_from_eul(eul);
  EXPECT_NEAR(quat_from_eul(0), true_quat(0), TOL);
  EXPECT_NEAR(quat_from_eul(1), true_quat(1), TOL);
  EXPECT_NEAR(quat_from_eul(2), true_quat(2), TOL);
  EXPECT_NEAR(quat_from_eul(3), true_quat(3), TOL);

  nt::RotationMatrix<double> dcm;
  dcm.setIdentity();
  nt::Quaternion<double> quat_from_dcm(dcm);
  EXPECT_NEAR(quat_from_dcm(0), 0.0, TOL);
  EXPECT_NEAR(quat_from_dcm(1), 0.0, TOL);
  EXPECT_NEAR(quat_from_dcm(2), 0.0, TOL);
  EXPECT_NEAR(quat_from_dcm(3), 1.0, TOL);
}

TEST(FunctionsTest, EulerAnglesFunctions) {
  // enforce euler angle creation with unnormalized angles
  nt::EulerAngles<double, nt::RotationOrder::YZX> eul;
  Eigen::Vector3d vec(nt::PI<double> + 0.1, -nt::PI<double> - 0.1, -nt::HALF_PI<double> - 0.1);
  eul = vec;
  eul.normalize();

  // get expected wrapped angles (YZX order)
  Eigen::Vector3d vec_norm;
  auto WrapAngle = [](double x, double y) {
    double half_y = 0.5 * y;
    return x - std::floor((x + half_y) / y) * y;
  };
  vec_norm(0) = WrapAngle(vec(0), nt::TWO_PI<double>);
  vec_norm(1) = WrapAngle(vec(1), nt::TWO_PI<double>);
  vec_norm(2) = WrapAngle(vec(2), nt::PI<double>);

  // test
  EXPECT_NEAR(eul(0), vec_norm(0), TOL);
  EXPECT_NEAR(eul(1), vec_norm(1), TOL);
  EXPECT_NEAR(eul(2), vec_norm(2), TOL);
}

TEST(FunctionsTest, DcmFunctions) {
  // Test RotationMatrix::Identity()
  nt::RotationMatrix<double> identity = nt::RotationMatrix<double>::Identity();
  ASSERT_TRUE(identity.isIdentity(TOL));

  // Test RotationMatrix::normalize()
  Eigen::Matrix3d mat;
  mat << 1.1, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 1.0;
  nt::RotationMatrix<double> dcm(mat);
  dcm.normalize();
  ASSERT_NEAR(dcm.determinant(), 1.0, TOL);

  // Test RotationMatrix::normalized()
  dcm = mat;
  nt::RotationMatrix<double> normalized_dcm = dcm.normalized();
  ASSERT_NEAR(normalized_dcm.determinant(), 1.0, TOL);
  ASSERT_FALSE(dcm.isIdentity(TOL));
  ASSERT_TRUE(normalized_dcm.isIdentity(TOL));
}

TEST(FunctionsTest, QuaternionFunctions) {
  // Test Quaternion::identity()
  nt::Quaternion<double> identity_quat = nt::Quaternion<double>::identity();
  ASSERT_NEAR(identity_quat.w(), 1.0, TOL);
  ASSERT_NEAR(identity_quat.x(), 0.0, TOL);
  ASSERT_NEAR(identity_quat.y(), 0.0, TOL);
  ASSERT_NEAR(identity_quat.z(), 0.0, TOL);

  // Test Quaternion::normalize() and normalized()
  nt::Quaternion<double> q_unnorm(1.0, 2.0, 3.0, 4.0);  // automatically normalized
  ASSERT_NEAR(q_unnorm.norm(), 1.0, TOL);
  nt::Quaternion<double> q_norm = q_unnorm.normalized();
  double norm1 = std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0 + 4.0 * 4.0);
  ASSERT_NEAR(q_norm.x(), 1.0 / norm1, TOL);
  ASSERT_NEAR(q_norm.y(), 2.0 / norm1, TOL);
  ASSERT_NEAR(q_norm.z(), 3.0 / norm1, TOL);
  ASSERT_NEAR(q_norm.w(), 4.0 / norm1, TOL);

  // Test Quaternion::conjugate() and conjugated()
  nt::Quaternion<double> q(0.1, 0.2, 0.3, 0.4);  // automatically normalized
  nt::Quaternion<double> q_conj = q.conjugated();
  double norm2 = std::sqrt(0.1 * 0.1 + 0.2 * 0.2 + 0.3 * 0.3 + 0.4 * 0.4);
  ASSERT_NEAR(q_conj.x(), -0.1 / norm2, TOL);
  ASSERT_NEAR(q_conj.y(), -0.2 / norm2, TOL);
  ASSERT_NEAR(q_conj.z(), -0.3 / norm2, TOL);
  ASSERT_NEAR(q_conj.w(), +0.4 / norm2, TOL);

  // Test Quaternion::inverse() and inversed()
  nt::Quaternion<double> unit_q(0.1, 0.2, 0.3, 0.4);  // automatically normalized
  nt::Quaternion<double> q_inv = unit_q.inversed();
  double norm3 = std::sqrt(0.1 * 0.1 + 0.2 * 0.2 + 0.3 * 0.3 + 0.4 * 0.4);
  ASSERT_NEAR(q_inv.x(), -0.1 / norm3, TOL);
  ASSERT_NEAR(q_inv.y(), -0.2 / norm3, TOL);
  ASSERT_NEAR(q_inv.z(), -0.3 / norm3, TOL);
  ASSERT_NEAR(q_inv.w(), +0.4 / norm3, TOL);
}

TEST(FunctionsTest, QuaternionProducts) {
  // Hamilton product (q * q)
  nt::Quaternion<double> q1(nt::SQRT_HALF<double>, nt::SQRT_HALF<double>, 0, 0);
  nt::Quaternion<double> q2(-nt::SQRT_HALF<double>, 0, 0, -nt::SQRT_HALF<double>);
  nt::Quaternion<double> q_result = q1 * q2;
  ASSERT_NEAR(q_result.x(), -0.5, TOL);
  ASSERT_NEAR(q_result.y(), -0.5, TOL);
  ASSERT_NEAR(q_result.z(), 0.5, TOL);
  ASSERT_NEAR(q_result.w(), 0.5, TOL);

  // Hamilton product (q *= q)
  nt::Quaternion<double> q1_inplace = q1;
  q1_inplace *= q2;
  ASSERT_NEAR(q_result.x(), -0.5, TOL);
  ASSERT_NEAR(q_result.y(), -0.5, TOL);
  ASSERT_NEAR(q_result.z(), 0.5, TOL);
  ASSERT_NEAR(q_result.w(), 0.5, TOL);

  // Rotation (q.rotate(v))
  nt::Quaternion<double> q_rot(nt::SQRT_HALF<double>, 0, 0, nt::SQRT_HALF<double>);
  Eigen::Vector3d vec(0, 1, 0);
  Eigen::Vector3d rotated_vec = q_rot.rotate(vec);
  ASSERT_NEAR(rotated_vec.x(), 0.0, TOL);
  ASSERT_NEAR(rotated_vec.y(), 0.0, TOL);
  ASSERT_NEAR(rotated_vec.z(), 1.0, TOL);
}