#include "root_functions.h"

using namespace std;


int main(int argc, char** argv) {

  plot_all_planes("X", "FWHM", true, argc, argv);

  //plot_all_planes_3D("Z","Y","FWHM",true, argc, argv);


  return 0;
}
