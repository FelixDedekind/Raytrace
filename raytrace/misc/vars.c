#include "vars.h"

// distances camera-lens, lens-object in meters
double d1 = 0.02, d2 = 0.05;

// FOV and resolution
// those can be made smaller LATER to save computation power while looking only at small angles
double theta_min = 0, phi_min = 0;
double theta_max = PI/200, phi_max = 2*PI;
// for efficiency reasons, theta_res should always be even.
// also, the variables should be chosen such that phi_res = theta_m/phi_m * theta_res
int theta_res = 50, phi_res = 1;