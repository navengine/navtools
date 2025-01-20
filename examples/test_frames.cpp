#include <iomanip>
#include <iostream>

#include "navtools/constants.hpp"
#include "navtools/frames.hpp"

const std::string RST = "\033[0m";  // reset
// const std::string RED = "\033[0;31m";   // red
// const std::string GRN = "\033[0;32m";   // green
// const std::string YEL = "\033[0;33m";   // yellow
// const std::string BLU = "\033[0;34m";   // blue
// const std::string MAG = "\033[0;35m";   // magenta
// const std::string CYN = "\033[0;36m";   // cyan
// const std::string WHT = "\033[0;37m";   // white
// const std::string BRED = "\033[1;31m";  // bold red
const std::string BGRN = "\033[1;32m";  // bold green
const std::string BYEL = "\033[1;33m";  // bold yellow
// const std::string BBLU = "\033[1;34m";  // bold blue
// const std::string BMAG = "\033[1;35m";  // bold magenta
const std::string BCYN = "\033[1;36m";  // bold cyan
// const std::string BWHT = "\033[1;37m";  // bold white

int main() {
  using namespace navtools;
  std::cout << std::setprecision(4) << std::endl;

  Eigen::Vector<double, 3> o = Eigen::Vector3d::Ones();
  std::cout << std::endl;

  // Starting with an ECEf position
  std::cout << BYEL << "#####* TESTING POSITION FRAME CONVERSIONS *#####" << RST << std::endl
            << std::endl;

  double dt = 1.0;
  Eigen::Vector<double, 3> starting_xyz{422596.629, -5362864.287, 3415493.797};
  Eigen::Vector<double, 3> moved_xyz = starting_xyz + o;
  Eigen::Vector<double, 3> eci =
      ecef2eci<double>(starting_xyz, dt);  //! {-2913450.7570, -4517890.5814, 3421252.4798}
  Eigen::Vector<double, 3> lla =
      ecef2lla<double>(starting_xyz);  //! {32.586279, -85.494372, 194.83}
  Eigen::Vector<double, 3> ned =
      ecef2ned<double>(moved_xyz, lla);  //! depends on correct lla {0,0,0}
  Eigen::Vector<double, 3> enu =
      ecef2enu<double>(moved_xyz, lla);  //! depends on correct lla {0,0,0}

  std::cout << BGRN << "Starting ECEF = " << starting_xyz.transpose() << RST << std::endl;
  std::cout << "ECI = " << eci.transpose() << std::endl;
  std::cout << "LLA = " << (lla.array() * LLA_RAD2DEG<double>.array()).transpose() << std::endl;
  std::cout << "NED = " << ned.transpose() << std::endl;
  std::cout << "ENU = " << enu.transpose() << std::endl << std::endl;

  // Rotate back to ECEF
  Eigen::Vector<double, 3> eciToxyz = eci2ecef<double>(eci, dt);
  Eigen::Vector<double, 3> llaToxyz = lla2ecef<double>(lla);
  Eigen::Vector<double, 3> nedToxyz = ned2ecef<double>(ned, lla);
  Eigen::Vector<double, 3> enuToxyz = enu2ecef<double>(enu, lla);

  std::cout << BGRN << "Starting ECEF = " << starting_xyz.transpose() << RST << std::endl;
  std::cout << "From ECI = " << eciToxyz.transpose() << std::endl;
  std::cout << "From LLA = " << llaToxyz.transpose() << std::endl;
  std::cout << "From NED = " << nedToxyz.transpose() << std::endl;
  std::cout << "From ENU = " << enuToxyz.transpose() << std::endl << std::endl;

  // Starting with NED Velocity
  std::cout << BYEL << "#####* TESTING VELOCITY FRAME CONVERSIONS *#####" << RST << std::endl
            << std::endl;

  Eigen::Vector<double, 3> starting_nedv{1.0, 2.0, 3.0};
  Eigen::Vector<double, 3> eciv = ned2eciv<double>(ned, starting_nedv, lla, dt);
  Eigen::Vector<double, 3> xyzv = ned2ecefv<double>(starting_nedv, lla);
  Eigen::Vector<double, 3> enuv = ned2enuDcm<double>() * starting_nedv;

  std::cout << BCYN << "Starting NEDV = " << starting_nedv.transpose() << RST << std::endl;
  std::cout << "ECIV = " << eciv.transpose() << std::endl;
  std::cout << "ECEFV = " << xyzv.transpose() << std::endl;
  std::cout << "ENUV = " << enuv.transpose() << std::endl << std::endl;

  Eigen::Vector<double, 3> eciTonedv = eci2nedv<double>(eci, eciv, lla, dt);
  Eigen::Vector<double, 3> xyzTonedv = ecef2nedv<double>(xyzv, lla);
  Eigen::Vector<double, 3> enuTonedv = enu2nedDcm<double>() * enuv;

  std::cout << BCYN << "Starting NEDV = " << starting_nedv.transpose() << RST << std::endl;
  std::cout << "From ECIV = " << eciTonedv.transpose() << std::endl;
  std::cout << "From ECEFV = " << xyzTonedv.transpose() << std::endl;
  std::cout << "From ENUV = " << enuTonedv.transpose() << std::endl << std::endl;

  std::cout << BYEL << "#####* DONE *#####" << RST << std::endl << std::endl;

  return 0;
}
