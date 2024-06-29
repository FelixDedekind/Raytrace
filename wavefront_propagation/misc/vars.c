#include "vars.h"

// distances camera-lens, lens-object in meters
double d1 = 0.05, d2 = 0.05;

double circle_radius = 0.001;
unsigned int circle_grid_num_phi = 1;
unsigned int circle_grid_num_r = 100;




double nAir = 1.00; // index of refraction of medium 1 (air)
double nWater =  1.33; // index of refraction of medium 2 (water)

const double wavelength = 500e-9;


double dot_radius = wavelength / 2.8 * 3;

double lens_radius = wavelength / 2.8 * 12;
unsigned int lens_grid_num_phi = 150;
unsigned int lens_grid_num_r = 150;