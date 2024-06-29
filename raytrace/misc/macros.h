#ifndef macros_defined
#define macros_defined

#define PI 3.14159265358979323846

#define FIELD_WIDTH 50
#define MAX_LINE_LENGTH 80

#define INDEX_OF_REFRACTION1 1.00 // index of refraction of medium 1 (air)
#define INDEX_OF_REFRACTION2 1.33 // index of refraction of medium 2 (water)

typedef struct {
    double x, y, z;
} Vertex;

typedef struct {
    Vertex* v1;
    Vertex* v2;
    Vertex* v3;
} Surface;

typedef struct {
    double x, y, z;
    double vx, vy, vz;
} Beam;

typedef struct {
    unsigned int r,g,b;
} Colour;


#endif 