#include <stdio.h> 
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#include "macros.h"



// distances camera-lens, lens-object in meters
extern double d1,d2;

extern double circle_radius;
extern unsigned int circle_grid_num_phi;
extern unsigned int circle_grid_num_r;



extern double nAir;
extern double nWater;

extern const double wavelength;

extern double dot_radius;

extern double lens_radius;
extern unsigned int lens_grid_num_phi;
extern unsigned int lens_grid_num_r;
