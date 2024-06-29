#include "vars.h"

// distances camera-lens, lens-object in meters
double d1 = 0.02, d2 = 0.05;

// FOV and resolution
// those can be made smaller LATER to save computation power while looking only at small angles
double theta_min = 0.00001, phi_min = 0.0000;
double theta_max = 0.02, phi_max = 2*PI;
// for efficiency reasons, theta_res should always be even.
// also, the variables should be chosen such that phi_res = theta_m/phi_m * theta_res
int theta_res = 5000, phi_res = 1;

double lens_radius = 0.005;

double circle_radius = 0.001;

double line_dist = 0.0001;
double line_width = 0.0001;