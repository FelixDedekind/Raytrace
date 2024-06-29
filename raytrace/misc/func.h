#include "vars.h"


void move_beam(Beam *beam, double s);
bool intersect_triangle(Beam *beam, Surface *surface, Vertex *intersection);
void refract(Beam* beam, Vertex* normal);


Colour evaluate_beam_circle(Beam *beam, double R);


void read_vertices(Vertex** vertices, const char *filename, int *num_thetas, int *num_phis);
void assign_surfaces(Surface **surfaces, Vertex *vertices, unsigned int num_thetas, unsigned int num_phis);
void transform_mirror_droplet(Vertex **vertices, double distance, int num_thetas, int num_phis);
void calculate_surface_normals(Surface* surfaces, Vertex* normals, int num_surfaces, Vertex* new_center);


void normalize_velocity(Beam *beam);
void initialize_beam(Beam *beam, int theta, int phi);
void print_beam(Beam *beam);
void print_evaluation(double theta, double phi, Colour colour, const char *filename);
void print_surfaces_to_file(Surface *surfaces, unsigned int num_surfaces, const char *filename, const char *type);
void process_indicator(char* s, char* status);
void export_surface_normals(Surface* surfaces, Vertex* normals, int num_surfaces);


void cross_product(Vertex* v1, Vertex* v2, Vertex* result);
void normalize(Vertex* v);
double dot_product(Vertex* v1, Vertex* v2);