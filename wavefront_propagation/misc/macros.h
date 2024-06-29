#ifndef macros_defined
#define macros_defined

#define PI 3.14159265358979323846

#define FIELD_WIDTH 50
#define MAX_LINE_LENGTH 80
#define ARROW_LENGTH 50

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
    double length;
} Beam;

typedef struct {
    unsigned int r,g,b;
} Colour;


#endif 