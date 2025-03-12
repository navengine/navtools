
#include <iomanip>
#include <iostream>

#include "navtools/attitude.hpp"
#include "navtools/constants.hpp"

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
// const std::string BCYN = "\033[1;36m";  // bold cyan
// const std::string BWHT = "\033[1;37m";  // bold white

int main() {
  using namespace navtools;
  std::cout << std::setprecision(4) << std::endl;

  Eigen::Vector<double, 3> starting_e{0.0, 5.0, -30.0};
  starting_e *= DEG2RAD<double>;

  Eigen::Vector<double, 4> qenu = euler2quat<double>(starting_e, false);
  Eigen::Vector<double, 4> qned = euler2quat<double>(starting_e, true);
  Eigen::Matrix<double, 3, 3> Renu = euler2dcm<double>(starting_e, false);
  Eigen::Matrix<double, 3, 3> Rned = euler2dcm<double>(starting_e, true);

  std::cout << BYEL << "#####* TESTING ATTITUDE TRANSFORMATIONS *#####" << RST << std::endl
            << std::endl;

  std::cout << BGRN << "Starting EULER = " << starting_e.transpose() << RST << std::endl;
  std::cout << "q_enu = " << qenu.transpose() << std::endl;
  std::cout << "q_ned = " << qned.transpose() << std::endl;
  std::cout << "C_enu = " << std::endl << Renu.transpose() << std::endl;
  std::cout << "C_ned = " << std::endl << Rned.transpose() << std::endl << std::endl;

  Eigen::Matrix<double, 3, 3> q2Renu = quat2dcm<double>(qenu);
  // dcmnorm(q2Renu);
  Eigen::Matrix<double, 3, 3> q2Rned = quat2dcm<double>(qned);
  // dcmnorm(q2Rned);
  Eigen::Vector<double, 4> R2qenu = dcm2quat<double>(Renu);
  // quatnorm(R2qenu);
  Eigen::Vector<double, 4> R2qned = dcm2quat<double>(Rned);
  // quatnorm(R2qned);

  std::cout << "C_enu from q_enu = " << std::endl << q2Renu << std::endl;
  std::cout << "C_ned from q_ned = " << std::endl << q2Rned << std::endl;
  std::cout << "q_enu from C_enu = " << R2qenu << std::endl;
  std::cout << "q_ned from C_ned = " << R2qned << std::endl << std::endl;

  Eigen::Vector<double, 3> e_Renu = dcm2euler<double>(q2Renu, false);
  Eigen::Vector<double, 3> e_Rned = dcm2euler<double>(q2Rned, true);
  Eigen::Vector<double, 3> e_qenu = quat2euler<double>(R2qenu, false);
  Eigen::Vector<double, 3> e_qned = quat2euler<double>(R2qned, true);

  std::cout << "e from q_enu = " << e_qenu.transpose() << std::endl;
  std::cout << "e from q_ned = " << e_qned.transpose() << std::endl;
  std::cout << "e from C_enu = " << e_Renu.transpose() << std::endl;
  std::cout << "e from C_ned = " << e_Rned.transpose() << std::endl << std::endl;

  return 0;
}