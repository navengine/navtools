
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

    Vec3<double> starting_e{0.0, 5.0, -30.0};
    starting_e *= DEG2RAD<double>;

    Vec4<double> qenu = euler2quat(starting_e, "enu");
    Vec4<double> qned = euler2quat(starting_e, "ned");
    Mat3x3<double> Renu = euler2dcm(starting_e, "enu");
    Mat3x3<double> Rned = euler2dcm(starting_e, "ned");

    std::cout << BYEL << "#####* TESTING ATTITUDE TRANSFORMATIONS *#####" << RST << std::endl
              << std::endl;

    std::cout << BGRN << "Starting EULER = " << starting_e.transpose() << RST << std::endl;
    std::cout << "q_enu = " << qenu.transpose() << std::endl;
    std::cout << "q_ned = " << qned.transpose() << std::endl;
    std::cout << "C_enu = " << std::endl << Renu.transpose() << std::endl;
    std::cout << "C_ned = " << std::endl << Rned.transpose() << std::endl << std::endl;

    Mat3x3<double> q2Renu = quat2dcm(qenu);
    // dcmnorm(q2Renu);
    Mat3x3<double> q2Rned = quat2dcm(qned);
    // dcmnorm(q2Rned);
    Vec4<double> R2qenu = dcm2quat(Renu);
    // quatnorm(R2qenu);
    Vec4<double> R2qned = dcm2quat(Rned);
    // quatnorm(R2qned);

    std::cout << "C_enu from q_enu = " << std::endl << q2Renu << std::endl;
    std::cout << "C_ned from q_ned = " << std::endl << q2Rned << std::endl;
    std::cout << "q_enu from C_enu = " << R2qenu << std::endl;
    std::cout << "q_ned from C_ned = " << R2qned << std::endl << std::endl;

    Vec3<double> e_Renu = dcm2euler(q2Renu, "enu");
    Vec3<double> e_Rned = dcm2euler(q2Rned, "ned");
    Vec3<double> e_qenu = quat2euler(R2qenu, "enu");
    Vec3<double> e_qned = quat2euler(R2qned, "ned");

    std::cout << "e from q_enu = " << e_qenu.transpose() << std::endl;
    std::cout << "e from q_ned = " << e_qned.transpose() << std::endl;
    std::cout << "e from C_enu = " << e_Renu.transpose() << std::endl;
    std::cout << "e from C_ned = " << e_Rned.transpose() << std::endl << std::endl;

    return 0;
}