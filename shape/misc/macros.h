#ifndef macros_defined
#define macros_defined

#define PI 3.14159265358979323846
#define g 9.81

#define FIELD_WIDTH 50


typedef struct {
    double x, y, z;
} Vertex;

typedef struct {
    Vertex* v1;
    Vertex* v2;
    Vertex* v3;
} Surface;

typedef struct {
    Surface *surfaces;
    Vertex **ptr_vertices;  //points to array of vertices! array saved externally. The array can be accessed via *(droplet->ptr_vertices)[0][ii]
    unsigned int num_surfaces;
    unsigned int num_vertices;
    unsigned int num_thetas;
    unsigned int num_phis;
    double density;
    double volume0;
    double sigma_to_air;
    double sigma_to_surface;
} Droplet;



#endif 