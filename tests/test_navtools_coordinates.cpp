#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <navtools/coordinates>

using EcefPos = nt::EcefCoord<nt::CoordType::POS, double>;
using EcefVel = nt::EcefCoord<nt::CoordType::VEL, double>;
using EcefAccel = nt::EcefCoord<nt::CoordType::ACCEL, double>;
using EcefAngVel = nt::EcefCoord<nt::CoordType::ANGVEL, double>;
using EciPos = nt::EciCoord<nt::CoordType::POS, double>;
using EciVel = nt::EciCoord<nt::CoordType::VEL, double>;
using EciAccel = nt::EciCoord<nt::CoordType::ACCEL, double>;
using EciAngVel = nt::EciCoord<nt::CoordType::ANGVEL, double>;
using NedPos = nt::TangentCoord<nt::CoordType::POS, nt::TangentFrame::NED, double>;
using NedVel = nt::TangentCoord<nt::CoordType::VEL, nt::TangentFrame::NED, double>;
using NedAccel = nt::TangentCoord<nt::CoordType::ACCEL, nt::TangentFrame::NED, double>;
using NedAngVel = nt::TangentCoord<nt::CoordType::ANGVEL, nt::TangentFrame::NED, double>;
using EnuPos = nt::TangentCoord<nt::CoordType::POS, nt::TangentFrame::ENU, double>;
using EnuVel = nt::TangentCoord<nt::CoordType::VEL, nt::TangentFrame::ENU, double>;
using EnuAccel = nt::TangentCoord<nt::CoordType::ACCEL, nt::TangentFrame::ENU, double>;
using EnuAngVel = nt::TangentCoord<nt::CoordType::ANGVEL, nt::TangentFrame::ENU, double>;
using LlaPosRad = nt::GeodeticCoord<false, double>;
using LlaPosDeg = nt::GeodeticCoord<true, double>;

constexpr double TOL = 1e-9;
const double dt = 3600.0;  // 1 hour

const LlaPosRad lla0(0.594, -2.063, 0.0);
const LlaPosRad lla_init(0.594314087066431, -2.06376826264174, 70.745);

const EcefPos r_ecef_init = nt::lla2ecef(lla_init);
const EciPos r_eci_init = nt::ecef2eci(dt, r_ecef_init);
const NedPos r_ned_init = nt::ecef2ned(r_ecef_init, lla0);
const EnuPos r_enu_init = nt::ecef2enu(r_ecef_init, lla0);

const EcefVel v_ecef_init(1.23, -4.56, 7.89);
const EciVel v_eci_init = nt::ecef2eciv(dt, v_ecef_init, r_ecef_init);
const NedVel v_ned_init = nt::ecef2nedv(v_ecef_init, lla0);
const EnuVel v_enu_init = nt::ecef2enuv(v_ecef_init, lla0);

const EcefAccel a_ecef_init(0.123, -0.456, 0.789);
const EciAccel a_eci_init = nt::ecef2ecia(dt, a_ecef_init, v_ecef_init, r_ecef_init);
const NedAccel a_ned_init = nt::ecef2neda(a_ecef_init, lla0);
const EnuAccel a_enu_init = nt::ecef2enua(a_ecef_init, lla0);

const EcefAngVel w_ecef_init(1.0, -0.2, 0.03);
const EciAngVel w_eci_init = nt::ecef2eciw(dt, w_ecef_init);
const NedAngVel w_ned_init = nt::ecef2nedw(w_ecef_init, lla0);
const EnuAngVel w_enu_init = nt::ecef2enuw(w_ecef_init, lla0);

TEST(CoordTransformations, PositionTransforms) {
  // LLA ->
  EcefPos lla_to_ecef = lla_init.to_ecef();
  EXPECT_NEAR(lla_to_ecef(0), r_ecef_init(0), TOL);
  EXPECT_NEAR(lla_to_ecef(1), r_ecef_init(1), TOL);
  EXPECT_NEAR(lla_to_ecef(2), r_ecef_init(2), TOL);
  EciPos lla_to_eci = lla_init.to_eci(dt);
  EXPECT_NEAR(lla_to_eci(0), r_eci_init(0), TOL);
  EXPECT_NEAR(lla_to_eci(1), r_eci_init(1), TOL);
  EXPECT_NEAR(lla_to_eci(2), r_eci_init(2), TOL);
  NedPos lla_to_ned = lla_init.to_tangent<nt::TangentFrame::NED>(lla0);
  EXPECT_NEAR(lla_to_ned(0), r_ned_init(0), TOL);
  EXPECT_NEAR(lla_to_ned(1), r_ned_init(1), TOL);
  EXPECT_NEAR(lla_to_ned(2), r_ned_init(2), TOL);
  EnuPos lla_to_enu = lla_init.to_tangent<nt::TangentFrame::ENU>(lla0);
  EXPECT_NEAR(lla_to_enu(0), r_enu_init(0), TOL);
  EXPECT_NEAR(lla_to_enu(1), r_enu_init(1), TOL);
  EXPECT_NEAR(lla_to_enu(2), r_enu_init(2), TOL);

  // LLA <-
  LlaPosRad lla_from_ecef = LlaPosRad::from_ecef(r_ecef_init);
  EXPECT_NEAR(lla_from_ecef(0), lla_init(0), TOL);
  EXPECT_NEAR(lla_from_ecef(1), lla_init(1), TOL);
  EXPECT_NEAR(lla_from_ecef(2), lla_init(2), TOL);
  LlaPosRad lla_from_eci = LlaPosRad::from_eci(dt, r_eci_init);
  EXPECT_NEAR(lla_from_eci(0), lla_init(0), TOL);
  EXPECT_NEAR(lla_from_eci(1), lla_init(1), TOL);
  EXPECT_NEAR(lla_from_eci(2), lla_init(2), TOL);
  LlaPosRad lla_from_ned = LlaPosRad::from_tangent<nt::TangentFrame::NED>(r_ned_init, lla0);
  EXPECT_NEAR(lla_from_ned(0), lla_init(0), TOL);
  EXPECT_NEAR(lla_from_ned(1), lla_init(1), TOL);
  EXPECT_NEAR(lla_from_ned(2), lla_init(2), TOL);
  LlaPosRad lla_from_enu = LlaPosRad::from_tangent<nt::TangentFrame::ENU>(r_enu_init, lla0);
  EXPECT_NEAR(lla_from_enu(0), lla_init(0), TOL);
  EXPECT_NEAR(lla_from_enu(1), lla_init(1), TOL);
  EXPECT_NEAR(lla_from_enu(2), lla_init(2), TOL);

  // ECEF ->
  LlaPosRad ecef_to_lla = r_ecef_init.to_lla();
  EXPECT_NEAR(ecef_to_lla(0), lla_init(0), TOL);
  EXPECT_NEAR(ecef_to_lla(1), lla_init(1), TOL);
  EXPECT_NEAR(ecef_to_lla(2), lla_init(2), TOL);
  EciPos ecef_to_eci = r_ecef_init.to_eci_pos(dt);
  EXPECT_NEAR(ecef_to_eci(0), r_eci_init(0), TOL);
  EXPECT_NEAR(ecef_to_eci(1), r_eci_init(1), TOL);
  EXPECT_NEAR(ecef_to_eci(2), r_eci_init(2), TOL);
  NedPos ecef_to_ned = r_ecef_init.to_tangent_pos<nt::TangentFrame::NED>(lla0);
  EXPECT_NEAR(ecef_to_ned(0), r_ned_init(0), TOL);
  EXPECT_NEAR(ecef_to_ned(1), r_ned_init(1), TOL);
  EXPECT_NEAR(ecef_to_ned(2), r_ned_init(2), TOL);
  EnuPos ecef_to_enu = r_ecef_init.to_tangent_pos<nt::TangentFrame::ENU>(lla0);
  EXPECT_NEAR(ecef_to_enu(0), r_enu_init(0), TOL);
  EXPECT_NEAR(ecef_to_enu(1), r_enu_init(1), TOL);
  EXPECT_NEAR(ecef_to_enu(2), r_enu_init(2), TOL);

  // ECEF <-
  EcefPos ecef_from_lla = EcefPos::from_lla(lla_init);
  EXPECT_NEAR(ecef_from_lla(0), r_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_lla(1), r_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_lla(2), r_ecef_init(2), TOL);
  EcefPos ecef_from_eci = EcefPos::from_eci_pos(dt, r_eci_init);
  EXPECT_NEAR(ecef_from_eci(0), r_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_eci(1), r_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_eci(2), r_ecef_init(2), TOL);
  EcefPos ecef_from_ned = EcefPos::from_tangent_pos<nt::TangentFrame::NED>(r_ned_init, lla0);
  EXPECT_NEAR(ecef_from_ned(0), r_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_ned(1), r_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_ned(2), r_ecef_init(2), TOL);
  EcefPos ecef_from_enu = EcefPos::from_tangent_pos<nt::TangentFrame::ENU>(r_enu_init, lla0);
  EXPECT_NEAR(ecef_from_enu(0), r_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_enu(1), r_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_enu(2), r_ecef_init(2), TOL);

  // ECI ->
  LlaPosRad eci_to_lla = r_eci_init.to_lla(dt);
  EXPECT_NEAR(eci_to_lla(0), lla_init(0), TOL);
  EXPECT_NEAR(eci_to_lla(1), lla_init(1), TOL);
  EXPECT_NEAR(eci_to_lla(2), lla_init(2), TOL);
  EcefPos eci_to_ecef = r_eci_init.to_ecef_pos(dt);
  EXPECT_NEAR(eci_to_ecef(0), r_ecef_init(0), TOL);
  EXPECT_NEAR(eci_to_ecef(1), r_ecef_init(1), TOL);
  EXPECT_NEAR(eci_to_ecef(2), r_ecef_init(2), TOL);
  NedPos eci_to_ned = r_eci_init.to_tangent_pos<nt::TangentFrame::NED>(dt, lla0);
  EXPECT_NEAR(eci_to_ned(0), r_ned_init(0), TOL);
  EXPECT_NEAR(eci_to_ned(1), r_ned_init(1), TOL);
  EXPECT_NEAR(eci_to_ned(2), r_ned_init(2), TOL);
  EnuPos eci_to_enu = r_eci_init.to_tangent_pos<nt::TangentFrame::ENU>(dt, lla0);
  EXPECT_NEAR(eci_to_enu(0), r_enu_init(0), TOL);
  EXPECT_NEAR(eci_to_enu(1), r_enu_init(1), TOL);
  EXPECT_NEAR(eci_to_enu(2), r_enu_init(2), TOL);

  // ECEF <-
  EciPos eci_from_lla = EciPos::from_lla(dt, lla_init);
  EXPECT_NEAR(eci_from_lla(0), r_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_lla(1), r_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_lla(2), r_eci_init(2), TOL);
  EciPos eci_from_eci = EciPos::from_ecef_pos(dt, r_ecef_init);
  EXPECT_NEAR(eci_from_eci(0), r_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_eci(1), r_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_eci(2), r_eci_init(2), TOL);
  EciPos eci_from_ned = EciPos::from_tangent_pos<nt::TangentFrame::NED>(dt, r_ned_init, lla0);
  EXPECT_NEAR(eci_from_ned(0), r_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_ned(1), r_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_ned(2), r_eci_init(2), TOL);
  EciPos eci_from_enu = EciPos::from_tangent_pos<nt::TangentFrame::ENU>(dt, r_enu_init, lla0);
  EXPECT_NEAR(eci_from_enu(0), r_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_enu(1), r_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_enu(2), r_eci_init(2), TOL);

  // NED ->
  LlaPosRad ned_to_lla = r_ned_init.to_lla(lla0);
  EXPECT_NEAR(ned_to_lla(0), lla_init(0), TOL);
  EXPECT_NEAR(ned_to_lla(1), lla_init(1), TOL);
  EXPECT_NEAR(ned_to_lla(2), lla_init(2), TOL);
  EciPos ned_to_eci = r_ned_init.to_eci_pos(dt, lla0);
  EXPECT_NEAR(ned_to_eci(0), r_eci_init(0), TOL);
  EXPECT_NEAR(ned_to_eci(1), r_eci_init(1), TOL);
  EXPECT_NEAR(ned_to_eci(2), r_eci_init(2), TOL);
  EcefPos ned_to_ecef = r_ned_init.to_ecef_pos(lla0);
  EXPECT_NEAR(ned_to_ecef(0), r_ecef_init(0), TOL);
  EXPECT_NEAR(ned_to_ecef(1), r_ecef_init(1), TOL);
  EXPECT_NEAR(ned_to_ecef(2), r_ecef_init(2), TOL);

  // NED <-
  NedPos ned_from_lla = NedPos::from_lla(lla_init, lla0);
  EXPECT_NEAR(ned_from_lla(0), r_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_lla(1), r_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_lla(2), r_ned_init(2), TOL);
  NedPos ned_from_eci = NedPos::from_eci_pos(dt, r_eci_init, lla0);
  EXPECT_NEAR(ned_from_eci(0), r_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_eci(1), r_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_eci(2), r_ned_init(2), TOL);
  NedPos ned_from_ecef = NedPos::from_ecef_pos(r_ecef_init, lla0);
  EXPECT_NEAR(ned_from_ecef(0), r_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_ecef(1), r_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_ecef(2), r_ned_init(2), TOL);

  // ENU ->
  LlaPosRad enu_to_lla = r_enu_init.to_lla(lla0);
  EXPECT_NEAR(enu_to_lla(0), lla_init(0), TOL);
  EXPECT_NEAR(enu_to_lla(1), lla_init(1), TOL);
  EXPECT_NEAR(enu_to_lla(2), lla_init(2), TOL);
  EciPos enu_to_eci = r_enu_init.to_eci_pos(dt, lla0);
  EXPECT_NEAR(enu_to_eci(0), r_eci_init(0), TOL);
  EXPECT_NEAR(enu_to_eci(1), r_eci_init(1), TOL);
  EXPECT_NEAR(enu_to_eci(2), r_eci_init(2), TOL);
  EcefPos enu_to_ecef = r_enu_init.to_ecef_pos(lla0);
  EXPECT_NEAR(enu_to_ecef(0), r_ecef_init(0), TOL);
  EXPECT_NEAR(enu_to_ecef(1), r_ecef_init(1), TOL);
  EXPECT_NEAR(enu_to_ecef(2), r_ecef_init(2), TOL);

  // ENU <-
  EnuPos enu_from_lla = EnuPos::from_lla(lla_init, lla0);
  EXPECT_NEAR(enu_from_lla(0), r_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_lla(1), r_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_lla(2), r_enu_init(2), TOL);
  EnuPos enu_from_eci = EnuPos::from_eci_pos(dt, r_eci_init, lla0);
  EXPECT_NEAR(enu_from_eci(0), r_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_eci(1), r_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_eci(2), r_enu_init(2), TOL);
  EnuPos enu_from_ecef = EnuPos::from_ecef_pos(r_ecef_init, lla0);
  EXPECT_NEAR(enu_from_ecef(0), r_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_ecef(1), r_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_ecef(2), r_enu_init(2), TOL);
}

TEST(CoordTransformations, VelocityTransforms) {
  // ECEF ->
  EciVel ecef_to_eci = v_ecef_init.to_eci_vel(dt, r_ecef_init);
  EXPECT_NEAR(ecef_to_eci(0), v_eci_init(0), TOL);
  EXPECT_NEAR(ecef_to_eci(1), v_eci_init(1), TOL);
  EXPECT_NEAR(ecef_to_eci(2), v_eci_init(2), TOL);
  NedVel ecef_to_ned = v_ecef_init.to_tangent_vel<nt::TangentFrame::NED>(lla0);
  EXPECT_NEAR(ecef_to_ned(0), v_ned_init(0), TOL);
  EXPECT_NEAR(ecef_to_ned(1), v_ned_init(1), TOL);
  EXPECT_NEAR(ecef_to_ned(2), v_ned_init(2), TOL);
  EnuVel ecef_to_enu = v_ecef_init.to_tangent_vel<nt::TangentFrame::ENU>(lla0);
  EXPECT_NEAR(ecef_to_enu(0), v_enu_init(0), TOL);
  EXPECT_NEAR(ecef_to_enu(1), v_enu_init(1), TOL);
  EXPECT_NEAR(ecef_to_enu(2), v_enu_init(2), TOL);

  // ECEF <-
  EcefVel ecef_from_eci = EcefVel::from_eci_vel(dt, v_eci_init, r_eci_init);
  EXPECT_NEAR(ecef_from_eci(0), v_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_eci(1), v_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_eci(2), v_ecef_init(2), TOL);
  EcefVel ecef_from_ned = EcefVel::from_tangent_vel<nt::TangentFrame::NED>(v_ned_init, lla0);
  EXPECT_NEAR(ecef_from_ned(0), v_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_ned(1), v_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_ned(2), v_ecef_init(2), TOL);
  EcefVel ecef_from_enu = EcefVel::from_tangent_vel<nt::TangentFrame::ENU>(v_enu_init, lla0);
  EXPECT_NEAR(ecef_from_enu(0), v_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_enu(1), v_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_enu(2), v_ecef_init(2), TOL);

  // ECI ->
  EcefVel eci_to_ecef = v_eci_init.to_ecef_vel(dt, r_eci_init);
  EXPECT_NEAR(eci_to_ecef(0), v_ecef_init(0), TOL);
  EXPECT_NEAR(eci_to_ecef(1), v_ecef_init(1), TOL);
  EXPECT_NEAR(eci_to_ecef(2), v_ecef_init(2), TOL);
  NedVel eci_to_ned = v_eci_init.to_tangent_vel<nt::TangentFrame::NED>(dt, r_eci_init, lla0);
  EXPECT_NEAR(eci_to_ned(0), v_ned_init(0), TOL);
  EXPECT_NEAR(eci_to_ned(1), v_ned_init(1), TOL);
  EXPECT_NEAR(eci_to_ned(2), v_ned_init(2), TOL);
  EnuVel eci_to_enu = v_eci_init.to_tangent_vel<nt::TangentFrame::ENU>(dt, r_eci_init, lla0);
  EXPECT_NEAR(eci_to_enu(0), v_enu_init(0), TOL);
  EXPECT_NEAR(eci_to_enu(1), v_enu_init(1), TOL);
  EXPECT_NEAR(eci_to_enu(2), v_enu_init(2), TOL);

  // ECI <-
  EciVel eci_from_ecef = EciVel::from_ecef_vel(dt, v_ecef_init, r_ecef_init);
  EXPECT_NEAR(eci_from_ecef(0), v_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_ecef(1), v_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_ecef(2), v_eci_init(2), TOL);
  EciVel eci_from_ned =
      EciVel::from_tangent_vel<nt::TangentFrame::NED>(dt, v_ned_init, r_ned_init, lla0);
  EXPECT_NEAR(eci_from_ned(0), v_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_ned(1), v_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_ned(2), v_eci_init(2), TOL);
  EciVel eci_from_enu =
      EciVel::from_tangent_vel<nt::TangentFrame::ENU>(dt, v_enu_init, r_enu_init, lla0);
  EXPECT_NEAR(eci_from_enu(0), v_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_enu(1), v_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_enu(2), v_eci_init(2), TOL);

  // NED ->
  EciVel ned_to_eci = v_ned_init.to_eci_vel(dt, r_ned_init, lla0);
  EXPECT_NEAR(ned_to_eci(0), v_eci_init(0), TOL);
  EXPECT_NEAR(ned_to_eci(1), v_eci_init(1), TOL);
  EXPECT_NEAR(ned_to_eci(2), v_eci_init(2), TOL);
  EcefVel ned_to_ecef = v_ned_init.to_ecef_vel(lla0);
  EXPECT_NEAR(ned_to_ecef(0), v_ecef_init(0), TOL);
  EXPECT_NEAR(ned_to_ecef(1), v_ecef_init(1), TOL);
  EXPECT_NEAR(ned_to_ecef(2), v_ecef_init(2), TOL);

  // NED <-
  NedVel ned_from_ecef = NedVel::from_ecef_vel(v_ecef_init, lla0);
  EXPECT_NEAR(ned_from_ecef(0), v_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_ecef(1), v_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_ecef(2), v_ned_init(2), TOL);
  NedVel ned_from_eci = NedVel::from_eci_vel(dt, v_eci_init, r_eci_init, lla0);
  EXPECT_NEAR(ned_from_eci(0), v_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_eci(1), v_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_eci(2), v_ned_init(2), TOL);

  // ENU ->
  EciVel enu_to_eci = v_enu_init.to_eci_vel(dt, r_enu_init, lla0);
  EXPECT_NEAR(enu_to_eci(0), v_eci_init(0), TOL);
  EXPECT_NEAR(enu_to_eci(1), v_eci_init(1), TOL);
  EXPECT_NEAR(enu_to_eci(2), v_eci_init(2), TOL);
  EcefVel enu_to_ecef = v_enu_init.to_ecef_vel(lla0);
  EXPECT_NEAR(enu_to_ecef(0), v_ecef_init(0), TOL);
  EXPECT_NEAR(enu_to_ecef(1), v_ecef_init(1), TOL);
  EXPECT_NEAR(enu_to_ecef(2), v_ecef_init(2), TOL);

  // ENU <-
  EnuVel enu_from_ecef = EnuVel::from_ecef_vel(v_ecef_init, lla0);
  EXPECT_NEAR(enu_from_ecef(0), v_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_ecef(1), v_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_ecef(2), v_enu_init(2), TOL);
  EnuVel enu_from_eci = EnuVel::from_eci_vel(dt, v_eci_init, r_eci_init, lla0);
  EXPECT_NEAR(enu_from_eci(0), v_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_eci(1), v_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_eci(2), v_enu_init(2), TOL);
}

TEST(CoordTransformations, AccelerationTransforms) {
  // ECEF ->
  EciAccel ecef_to_eci = a_ecef_init.to_eci_accel(dt, v_ecef_init, r_ecef_init);
  EXPECT_NEAR(ecef_to_eci(0), a_eci_init(0), TOL);
  EXPECT_NEAR(ecef_to_eci(1), a_eci_init(1), TOL);
  EXPECT_NEAR(ecef_to_eci(2), a_eci_init(2), TOL);
  NedAccel ecef_to_ned = a_ecef_init.to_tangent_accel<nt::TangentFrame::NED>(lla0);
  EXPECT_NEAR(ecef_to_ned(0), a_ned_init(0), TOL);
  EXPECT_NEAR(ecef_to_ned(1), a_ned_init(1), TOL);
  EXPECT_NEAR(ecef_to_ned(2), a_ned_init(2), TOL);
  EnuAccel ecef_to_enu = a_ecef_init.to_tangent_accel<nt::TangentFrame::ENU>(lla0);
  EXPECT_NEAR(ecef_to_enu(0), a_enu_init(0), TOL);
  EXPECT_NEAR(ecef_to_enu(1), a_enu_init(1), TOL);
  EXPECT_NEAR(ecef_to_enu(2), a_enu_init(2), TOL);

  // ECEF <-
  EcefAccel ecef_from_eci = EcefAccel::from_eci_accel(dt, a_eci_init, v_eci_init, r_eci_init);
  EXPECT_NEAR(ecef_from_eci(0), a_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_eci(1), a_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_eci(2), a_ecef_init(2), TOL);
  EcefAccel ecef_from_ned = EcefAccel::from_tangent_accel<nt::TangentFrame::NED>(a_ned_init, lla0);
  EXPECT_NEAR(ecef_from_ned(0), a_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_ned(1), a_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_ned(2), a_ecef_init(2), TOL);
  EcefAccel ecef_from_enu = EcefAccel::from_tangent_accel<nt::TangentFrame::ENU>(a_enu_init, lla0);
  EXPECT_NEAR(ecef_from_enu(0), a_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_enu(1), a_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_enu(2), a_ecef_init(2), TOL);

  // ECI ->
  EcefAccel eci_to_ecef = a_eci_init.to_ecef_accel(dt, v_eci_init, r_eci_init);
  EXPECT_NEAR(eci_to_ecef(0), a_ecef_init(0), TOL);
  EXPECT_NEAR(eci_to_ecef(1), a_ecef_init(1), TOL);
  EXPECT_NEAR(eci_to_ecef(2), a_ecef_init(2), TOL);
  NedAccel eci_to_ned =
      a_eci_init.to_tangent_accel<nt::TangentFrame::NED>(dt, v_eci_init, r_eci_init, lla0);
  EXPECT_NEAR(eci_to_ned(0), a_ned_init(0), TOL);
  EXPECT_NEAR(eci_to_ned(1), a_ned_init(1), TOL);
  EXPECT_NEAR(eci_to_ned(2), a_ned_init(2), TOL);
  EnuAccel eci_to_enu =
      a_eci_init.to_tangent_accel<nt::TangentFrame::ENU>(dt, v_eci_init, r_eci_init, lla0);
  EXPECT_NEAR(eci_to_enu(0), a_enu_init(0), TOL);
  EXPECT_NEAR(eci_to_enu(1), a_enu_init(1), TOL);
  EXPECT_NEAR(eci_to_enu(2), a_enu_init(2), TOL);

  // ECI <-
  EciAccel eci_from_ecef = EciAccel::from_ecef_accel(dt, a_ecef_init, v_ecef_init, r_ecef_init);
  EXPECT_NEAR(eci_from_ecef(0), a_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_ecef(1), a_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_ecef(2), a_eci_init(2), TOL);
  EciAccel eci_from_ned = EciAccel::from_tangent_accel<nt::TangentFrame::NED>(
      dt, a_ned_init, v_ned_init, r_ned_init, lla0);
  EXPECT_NEAR(eci_from_ned(0), a_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_ned(1), a_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_ned(2), a_eci_init(2), TOL);
  EciAccel eci_from_enu = EciAccel::from_tangent_accel<nt::TangentFrame::ENU>(
      dt, a_enu_init, v_enu_init, r_enu_init, lla0);
  EXPECT_NEAR(eci_from_enu(0), a_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_enu(1), a_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_enu(2), a_eci_init(2), TOL);

  // NED ->
  EciAccel ned_to_eci = a_ned_init.to_eci_accel(dt, v_ned_init, r_ned_init, lla0);
  EXPECT_NEAR(ned_to_eci(0), a_eci_init(0), TOL);
  EXPECT_NEAR(ned_to_eci(1), a_eci_init(1), TOL);
  EXPECT_NEAR(ned_to_eci(2), a_eci_init(2), TOL);
  EcefAccel ned_to_ecef = a_ned_init.to_ecef_accel(lla0);
  EXPECT_NEAR(ned_to_ecef(0), a_ecef_init(0), TOL);
  EXPECT_NEAR(ned_to_ecef(1), a_ecef_init(1), TOL);
  EXPECT_NEAR(ned_to_ecef(2), a_ecef_init(2), TOL);

  // NED <-
  NedAccel ned_from_ecef = NedAccel::from_ecef_accel(a_ecef_init, lla0);
  EXPECT_NEAR(ned_from_ecef(0), a_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_ecef(1), a_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_ecef(2), a_ned_init(2), TOL);
  NedAccel ned_from_eci = NedAccel::from_eci_accel(dt, a_eci_init, v_eci_init, r_eci_init, lla0);
  EXPECT_NEAR(ned_from_eci(0), a_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_eci(1), a_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_eci(2), a_ned_init(2), TOL);

  // ENU ->
  EciAccel enu_to_eci = a_enu_init.to_eci_accel(dt, v_enu_init, r_enu_init, lla0);
  EXPECT_NEAR(enu_to_eci(0), a_eci_init(0), TOL);
  EXPECT_NEAR(enu_to_eci(1), a_eci_init(1), TOL);
  EXPECT_NEAR(enu_to_eci(2), a_eci_init(2), TOL);
  EcefAccel enu_to_ecef = a_enu_init.to_ecef_accel(lla0);
  EXPECT_NEAR(enu_to_ecef(0), a_ecef_init(0), TOL);
  EXPECT_NEAR(enu_to_ecef(1), a_ecef_init(1), TOL);
  EXPECT_NEAR(enu_to_ecef(2), a_ecef_init(2), TOL);

  // ENU <-
  EnuAccel enu_from_ecef = EnuAccel::from_ecef_accel(a_ecef_init, lla0);
  EXPECT_NEAR(enu_from_ecef(0), a_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_ecef(1), a_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_ecef(2), a_enu_init(2), TOL);
  EnuAccel enu_from_eci = EnuAccel::from_eci_accel(dt, a_eci_init, v_eci_init, r_eci_init, lla0);
  EXPECT_NEAR(enu_from_eci(0), a_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_eci(1), a_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_eci(2), a_enu_init(2), TOL);
}

TEST(CoordTransformations, AngularVelocityTransforms) {
  // ECEF ->
  EciAngVel ecef_to_eci = w_ecef_init.to_eci_angvel(dt);
  EXPECT_NEAR(ecef_to_eci(0), w_eci_init(0), TOL);
  EXPECT_NEAR(ecef_to_eci(1), w_eci_init(1), TOL);
  EXPECT_NEAR(ecef_to_eci(2), w_eci_init(2), TOL);
  NedAngVel ecef_to_ned = w_ecef_init.to_tangent_angvel<nt::TangentFrame::NED>(lla0);
  EXPECT_NEAR(ecef_to_ned(0), w_ned_init(0), TOL);
  EXPECT_NEAR(ecef_to_ned(1), w_ned_init(1), TOL);
  EXPECT_NEAR(ecef_to_ned(2), w_ned_init(2), TOL);
  EnuAngVel ecef_to_enu = w_ecef_init.to_tangent_angvel<nt::TangentFrame::ENU>(lla0);
  EXPECT_NEAR(ecef_to_enu(0), w_enu_init(0), TOL);
  EXPECT_NEAR(ecef_to_enu(1), w_enu_init(1), TOL);
  EXPECT_NEAR(ecef_to_enu(2), w_enu_init(2), TOL);

  // ECEF <-
  EcefAngVel ecef_from_eci = EcefAngVel::from_eci_angvel(dt, w_eci_init);
  EXPECT_NEAR(ecef_from_eci(0), w_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_eci(1), w_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_eci(2), w_ecef_init(2), TOL);
  EcefAngVel ecef_from_ned =
      EcefAngVel::from_tangent_angvel<nt::TangentFrame::NED>(w_ned_init, lla0);
  EXPECT_NEAR(ecef_from_ned(0), w_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_ned(1), w_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_ned(2), w_ecef_init(2), TOL);
  EcefAngVel ecef_from_enu =
      EcefAngVel::from_tangent_angvel<nt::TangentFrame::ENU>(w_enu_init, lla0);
  EXPECT_NEAR(ecef_from_enu(0), w_ecef_init(0), TOL);
  EXPECT_NEAR(ecef_from_enu(1), w_ecef_init(1), TOL);
  EXPECT_NEAR(ecef_from_enu(2), w_ecef_init(2), TOL);

  // ECI ->
  EcefAngVel eci_to_ecef = w_eci_init.to_ecef_angvel(dt);
  EXPECT_NEAR(eci_to_ecef(0), w_ecef_init(0), TOL);
  EXPECT_NEAR(eci_to_ecef(1), w_ecef_init(1), TOL);
  EXPECT_NEAR(eci_to_ecef(2), w_ecef_init(2), TOL);
  NedAngVel eci_to_ned = w_eci_init.to_tangent_angvel<nt::TangentFrame::NED>(dt, lla0);
  EXPECT_NEAR(eci_to_ned(0), w_ned_init(0), TOL);
  EXPECT_NEAR(eci_to_ned(1), w_ned_init(1), TOL);
  EXPECT_NEAR(eci_to_ned(2), w_ned_init(2), TOL);
  EnuAngVel eci_to_enu = w_eci_init.to_tangent_angvel<nt::TangentFrame::ENU>(dt, lla0);
  EXPECT_NEAR(eci_to_enu(0), w_enu_init(0), TOL);
  EXPECT_NEAR(eci_to_enu(1), w_enu_init(1), TOL);
  EXPECT_NEAR(eci_to_enu(2), w_enu_init(2), TOL);

  // ECI <-
  EciAngVel eci_from_ecef = EciAngVel::from_ecef_angvel(dt, w_ecef_init);
  EXPECT_NEAR(eci_from_ecef(0), w_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_ecef(1), w_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_ecef(2), w_eci_init(2), TOL);
  EciAngVel eci_from_ned =
      EciAngVel::from_tangent_angvel<nt::TangentFrame::NED>(dt, w_ned_init, lla0);
  EXPECT_NEAR(eci_from_ned(0), w_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_ned(1), w_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_ned(2), w_eci_init(2), TOL);
  EciAngVel eci_from_enu =
      EciAngVel::from_tangent_angvel<nt::TangentFrame::ENU>(dt, w_enu_init, lla0);
  EXPECT_NEAR(eci_from_enu(0), w_eci_init(0), TOL);
  EXPECT_NEAR(eci_from_enu(1), w_eci_init(1), TOL);
  EXPECT_NEAR(eci_from_enu(2), w_eci_init(2), TOL);

  // NED ->
  EciAngVel ned_to_eci = w_ned_init.to_eci_angvel(dt, lla0);
  EXPECT_NEAR(ned_to_eci(0), w_eci_init(0), TOL);
  EXPECT_NEAR(ned_to_eci(1), w_eci_init(1), TOL);
  EXPECT_NEAR(ned_to_eci(2), w_eci_init(2), TOL);
  EcefAngVel ned_to_ecef = w_ned_init.to_ecef_angvel(lla0);
  EXPECT_NEAR(ned_to_ecef(0), w_ecef_init(0), TOL);
  EXPECT_NEAR(ned_to_ecef(1), w_ecef_init(1), TOL);
  EXPECT_NEAR(ned_to_ecef(2), w_ecef_init(2), TOL);

  // NED <-
  NedAngVel ned_from_ecef = NedAngVel::from_ecef_angvel(w_ecef_init, lla0);
  EXPECT_NEAR(ned_from_ecef(0), w_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_ecef(1), w_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_ecef(2), w_ned_init(2), TOL);
  NedAngVel ned_from_eci = NedAngVel::from_eci_angvel(dt, w_eci_init, lla0);
  EXPECT_NEAR(ned_from_eci(0), w_ned_init(0), TOL);
  EXPECT_NEAR(ned_from_eci(1), w_ned_init(1), TOL);
  EXPECT_NEAR(ned_from_eci(2), w_ned_init(2), TOL);

  // ENU ->
  EciAngVel enu_to_eci = w_enu_init.to_eci_angvel(dt, lla0);
  EXPECT_NEAR(enu_to_eci(0), w_eci_init(0), TOL);
  EXPECT_NEAR(enu_to_eci(1), w_eci_init(1), TOL);
  EXPECT_NEAR(enu_to_eci(2), w_eci_init(2), TOL);
  EcefAngVel enu_to_ecef = w_enu_init.to_ecef_angvel(lla0);
  EXPECT_NEAR(enu_to_ecef(0), w_ecef_init(0), TOL);
  EXPECT_NEAR(enu_to_ecef(1), w_ecef_init(1), TOL);
  EXPECT_NEAR(enu_to_ecef(2), w_ecef_init(2), TOL);

  // ENU <-
  EnuAngVel enu_from_ecef = EnuAngVel::from_ecef_angvel(w_ecef_init, lla0);
  EXPECT_NEAR(enu_from_ecef(0), w_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_ecef(1), w_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_ecef(2), w_enu_init(2), TOL);
  EnuAngVel enu_from_eci = EnuAngVel::from_eci_angvel(dt, w_eci_init, lla0);
  EXPECT_NEAR(enu_from_eci(0), w_enu_init(0), TOL);
  EXPECT_NEAR(enu_from_eci(1), w_enu_init(1), TOL);
  EXPECT_NEAR(enu_from_eci(2), w_enu_init(2), TOL);
}