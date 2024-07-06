#include "func.h"
#include "vars.h"


// --------------------------------------------------- raytrace functions

void move_beam(Beam *beam, double s) {
    beam->x += beam->vx * s;
    beam->y += beam->vy * s;
    beam->z += beam->vz * s;
}


bool intersect_triangle(Beam *beam, Surface *surface, Vertex *intersection) {
    double scale_factor = 1000.0; // Adjust this value as needed

    // Create temporary scaled vertices
    Vertex v1_scaled = {surface->v1->x * scale_factor, surface->v1->y * scale_factor, surface->v1->z * scale_factor};
    Vertex v2_scaled = {surface->v2->x * scale_factor, surface->v2->y * scale_factor, surface->v2->z * scale_factor};
    Vertex v3_scaled = {surface->v3->x * scale_factor, surface->v3->y * scale_factor, surface->v3->z * scale_factor};

    // Create a temporary scaled beam
    Beam beam_scaled = {beam->x * scale_factor, beam->y * scale_factor, beam->z * scale_factor, beam->vx, beam->vy, beam->vz};

    Vertex edge1, edge2, h, s, q;
    double a,f,u,v;
    edge1.x = v2_scaled.x - v1_scaled.x;
    edge1.y = v2_scaled.y - v1_scaled.y;
    edge1.z = v2_scaled.z - v1_scaled.z;

    edge2.x = v3_scaled.x - v1_scaled.x;
    edge2.y = v3_scaled.y - v1_scaled.y;
    edge2.z = v3_scaled.z - v1_scaled.z;

    h.x = beam_scaled.vy * edge2.z - beam_scaled.vz * edge2.y;
    h.y = beam_scaled.vz * edge2.x - beam_scaled.vx * edge2.z;
    h.z = beam_scaled.vx * edge2.y - beam_scaled.vy * edge2.x;

    a = edge1.x * h.x + edge1.y * h.y + edge1.z * h.z;

    if (a > -0.00001 && a < 0.00001)
        return(false);

    f = 1/a;
    s.x = beam_scaled.x - v1_scaled.x;
    s.y = beam_scaled.y - v1_scaled.y;
    s.z = beam_scaled.z - v1_scaled.z;

    u = f * (s.x * h.x + s.y * h.y + s.z * h.z);

    if (u < 0.0 || u > 1.0)
        return(false);

    q.x = s.y * edge1.z - s.z * edge1.y;
    q.y = s.z * edge1.x - s.x * edge1.z;
    q.z = s.x * edge1.y - s.y * edge1.x;

    v = f * (beam_scaled.vx * q.x + beam_scaled.vy * q.y + beam_scaled.vz * q.z);

    if (v < 0.0 || u + v > 1.0)
        return(false);

    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * (edge2.x * q.x + edge2.y * q.y + edge2.z * q.z);

    if (t > 0.00001) { // ray intersection
        // Compute the intersection point and assign it to the intersection vertex
        intersection->x = beam_scaled.x + beam_scaled.vx * t;
        intersection->y = beam_scaled.y + beam_scaled.vy * t;
        intersection->z = beam_scaled.z + beam_scaled.vz * t;

        // Scale down the intersection point to the original scale
        intersection->x /= scale_factor;
        intersection->y /= scale_factor;
        intersection->z /= scale_factor;

        return(true);
    } else { // This means that there is a line intersection but not a ray intersection.
        return (false);
    }
}


void refract(Beam* beam, Vertex* normal) {
    double dot_product = beam->vx * normal->x + beam->vy * normal->y + beam->vz * normal->z;
    double k = 1.0 - INDEX_OF_REFRACTION1 / INDEX_OF_REFRACTION2 * (1.0 - dot_product * dot_product);

    if (k < 0.0) {
        // Total internal reflection
        beam->vx = 2.0 * dot_product * normal->x - beam->vx;
        beam->vy = 2.0 * dot_product * normal->y - beam->vy;
        beam->vz = 2.0 * dot_product * normal->z - beam->vz;
    } else {
        // Refraction
        double scale = INDEX_OF_REFRACTION1 / INDEX_OF_REFRACTION2;
        beam->vx = scale * beam->vx - (scale * dot_product + sqrt(k)) * normal->x;
        beam->vy = scale * beam->vy - (scale * dot_product + sqrt(k)) * normal->y;
        beam->vz = scale * beam->vz - (scale * dot_product + sqrt(k)) * normal->z;
    }
}







// --------------------------------------------------- evaluation functions

Colour evaluate_beam_circle(Beam *beam, double R) {
    Colour colour;
    colour.r = 0;
    colour.g = 0;
    colour.b = 0;
    double r = sqrt(beam->x*beam->x + beam->y*beam->y);
    if(r < R)   colour.r = 255;
    return colour;
}

Colour evaluate_beam_lines(Beam *beam) {
    Colour colour;
    colour.r = 0;
    colour.g = 0;
    colour.b = 0;
    
    if((-line_width - line_dist/2 < beam->x && beam->x < -line_dist/2) || (line_width + line_dist/2 > beam->x && beam->x > line_dist/2))
        colour.r = 255;

    return colour;
}



// --------------------------------------------------- data functions

// read in positions of vertices from file
// ORDER OF VERTICES IS IMPORTANT from now on, I don't want to define all surfaces here.
void read_vertices(Vertex** vertices, const char *filename, int *num_thetas, int *num_phis) {
    // open file
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Could not open file\n");
        return;
    }

    // scan for num_thetas and num_phi
    fscanf(file, "%d %d", num_thetas, num_phis);


    // allocate memory for vertices
    unsigned int num_vertices = (*num_phis * (*num_thetas - 1) + 1);
    *vertices = (Vertex*) malloc(num_vertices * sizeof(Vertex));
    if (*vertices == NULL) {
        printf("Memory allocation failed\n");
        return;
    }

    //helper variables
    double a, b, c;
    int index = 0;

    // read in rest of lines
    while (fscanf(file, "%lf %lf %lf", &a, &b, &c) == 3) {
        (*vertices)[index].x = a;
        (*vertices)[index].y = b;
        (*vertices)[index].z = c;
        index++;
    }

    fclose(file);
}


void read_vertices_2d(Vertex** vertices, const char *filename, int *num_thetas, int *num_phis) {
    // open file
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Could not open file\n");
        return;
    }

    // scan for num_thetas and num_phi
    fscanf(file, "%d %d", num_thetas, num_phis);


    // allocate memory for vertices
    unsigned int num_vertices = (*num_phis * (*num_thetas - 1) + 1);
    *vertices = (Vertex*) malloc(num_vertices * sizeof(Vertex));
    if (*vertices == NULL) {
        printf("Memory allocation failed\n");
        return;
    }

    //helper variables
    double a, b, c;
    int index = 0;

    // read in rest of lines
    while (fscanf(file, "%lf %lf %lf", &a, &b, &c) == 3) {
        if(index == 0) 
        {
            (*vertices)[index].x = a;
            (*vertices)[index].y = b;
            (*vertices)[index].z = c;
            index++;
        }
        else 
        {
            for(int ii = 0; ii < *num_phis; ii++)
            {
                (*vertices)[index].x = a*cos((double)ii/ *num_phis*2*PI);
                (*vertices)[index].y = a*sin((double)ii/ *num_phis*2*PI);
                (*vertices)[index].z = c;
                index++;
            }
        }
        
    }

    fclose(file);
}


void assign_surfaces(Surface **surfaces, Vertex *vertices, unsigned int num_thetas, unsigned int num_phis) {
    unsigned int num_surfaces = (num_thetas - 2) * 2 * num_phis + num_phis;
    *surfaces = malloc(num_surfaces * sizeof(Surface));
    int tt, pp;

    for (pp = 0; pp < num_phis; pp++) {
        // Create one surface for each grid cell
        (*surfaces)[pp].v1 = &(vertices)[0];
        (*surfaces)[pp].v2 = &(vertices)[(pp + 1)];
        (*surfaces)[pp].v3 = &(vertices)[(pp + 1) % num_phis + 1];
    }

    unsigned int vertex_offset = 1;
    unsigned int surface_offset = num_phis;

    // then the middle surfaces
    for (tt = 0; tt < num_thetas - 2; tt++) {
        for (pp = 0; pp < num_phis; pp++) {
            // printf("surface %i \n", 2 * tt * num_phis + surface_offset + 2 * pp +1);
            // printf("vertex: %i \n", vertex_offset + (tt+1)*num_phis + ((pp + 1) % num_phis));
            // Create two surfaces for each grid cell
            (*surfaces)[2 * tt * num_phis + surface_offset + 2 * pp].v1 = &(vertices)[vertex_offset + tt*num_phis + pp];
            (*surfaces)[2 * tt * num_phis + surface_offset + 2 * pp].v2 = &(vertices)[vertex_offset + tt*num_phis + ((pp + 1) % num_phis)];
            (*surfaces)[2 * tt * num_phis + surface_offset + 2 * pp].v3 = &(vertices)[vertex_offset + (tt+1)*num_phis + pp];

            (*surfaces)[2 * tt * num_phis + surface_offset + 2 * pp + 1].v1 = &(vertices)[vertex_offset + tt*num_phis + ((pp + 1) % num_phis)];
            (*surfaces)[2 * tt * num_phis + surface_offset + 2 * pp + 1].v2 = &(vertices)[vertex_offset + (tt+1)*num_phis + pp];
            (*surfaces)[2 * tt * num_phis + surface_offset + 2 * pp + 1].v3 = &(vertices)[vertex_offset + (tt+1)*num_phis + ((pp + 1) % num_phis)];
        }
    }
}


void transform_mirror_droplet(Vertex **vertices, double distance, int num_thetas, int num_phis) {
    unsigned int num_vertices = (num_phis * (num_thetas - 1) + 1);  // this only works for half-sphere
    for(int ii = 0; ii < num_vertices; ii++) {
        (*vertices)[ii].z = distance - (*vertices)[ii].z;
    }
}


void calculate_surface_normals(Surface* surfaces, Vertex* normals, int num_surfaces, Vertex* new_center) {
    for (int i = 0; i < num_surfaces; i++) {
        // Calculate the normal vector
        Vertex v1 = {surfaces[i].v2->x - surfaces[i].v1->x, surfaces[i].v2->y - surfaces[i].v1->y, surfaces[i].v2->z - surfaces[i].v1->z};
        Vertex v2 = {surfaces[i].v3->x - surfaces[i].v1->x, surfaces[i].v3->y - surfaces[i].v1->y, surfaces[i].v3->z - surfaces[i].v1->z};
        cross_product(&v1, &v2, &normals[i]);
        normalize(&normals[i]);

        // Ensure the normal points outwards
        Vertex radial = {surfaces[i].v1->x - new_center->x, surfaces[i].v1->y - new_center->y, surfaces[i].v1->z - new_center->z};  // vector from new center to a point on the surface
        if (dot_product(&normals[i], &radial) < 0) {
            normals[i].x = -normals[i].x;
            normals[i].y = -normals[i].y;
            normals[i].z = -normals[i].z;
        }
    }
}






// --------------------------------------------------- helper functions

void normalize_velocity(Beam *beam) {
    double magnitude = sqrt(beam->vx * beam->vx + beam->vy * beam->vy + beam->vz * beam->vz);
    if (magnitude != 0) { // to avoid division by zero
        beam->vx /= magnitude;
        beam->vy /= magnitude;
        beam->vz /= magnitude;
    }
}


// want to sample angles uniformly
void initialize_beam(Beam *beam, int theta, int phi) {
    // Initialize the position of the beam at the camera (0,0,0)
    beam->x = 0;
    beam->y = 0;
    beam->z = 0;

    // Calculate the step size for theta and phi
    double theta_step = (theta_max - theta_min) / (double)theta_res;
    double phi_step = (phi_max - phi_min) / (double)phi_res;

    // Calculate the actual theta and phi values
    double theta_val = theta_min + (double)theta * theta_step;
    double phi_val = phi_min + ((double)phi + 0.5) * phi_step; // 0.5 is so that it hits centre for phi_res = 1

    // Convert spherical coordinates to Cartesian coordinates
    beam->vx = sin(theta_val) * cos(phi_val);
    beam->vy = sin(theta_val) * sin(phi_val);
    beam->vz = cos(theta_val);
}


void print_beam(Beam *beam) {
    printf("beam position: \t %f %f %f \nbeam direction:\t %f %f %f \n", beam->x, beam->y, beam->z, beam->vx, beam->vy, beam->vz);
}


void print_evaluation(double theta, double phi, Colour colour, const char *filename) {
    FILE *file = fopen(filename, "a");
    fprintf(file, "%0.20f %0.20f %i %i %i \n", theta, phi, colour.r, colour.g, colour.b);
    fclose(file);
}


void print_surfaces_to_file(Surface *surfaces, unsigned int num_surfaces, const char *filename, const char *type) {
    FILE *file = fopen(filename, type);
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    // printf("num surfaces: %i \n", droplet->num_surfaces);

    for (int i = 0; i < num_surfaces; i++) {
        // printf("trying to write to file: surface number %i \n", i);
        Surface *s = &surfaces[i];
        fprintf(file, "%0.10f %0.10f %0.10f , %0.10f %0.10f %0.10f , %0.10f %0.10f %0.10f \n",
                s->v1->x, s->v1->y, s->v1->z,
                s->v2->x, s->v2->y, s->v2->z,
                s->v3->x, s->v3->y, s->v3->z);
    }

    fclose(file);
}


void process_indicator(char* s, char* status){
    printf("%-*s%s\n", FIELD_WIDTH, s ,status);
}


void export_surface_normals(Surface* surfaces, Vertex* normals, int num_surfaces) {
    FILE* file = fopen("surface_normals.txt", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < num_surfaces; i++) {
        // Calculate the center of the surface
        Vertex center = {(surfaces[i].v1->x + surfaces[i].v2->x + surfaces[i].v3->x) / 3.0,
                         (surfaces[i].v1->y + surfaces[i].v2->y + surfaces[i].v3->y) / 3.0,
                         (surfaces[i].v1->z + surfaces[i].v2->z + surfaces[i].v3->z) / 3.0};

        // Write the center and normal to the file
        fprintf(file, "%f %f %f, %f %f %f\n", center.x, center.y, center.z, normals[i].x, normals[i].y, normals[i].z);
    }

    fclose(file);
}











// --------------------------------------------------- vertex functions

// Function to calculate the cross product of two vectors
void cross_product(Vertex* v1, Vertex* v2, Vertex* result) {
    result->x = v1->y * v2->z - v1->z * v2->y;
    result->y = v1->z * v2->x - v1->x * v2->z;
    result->z = v1->x * v2->y - v1->y * v2->x;
}

// Function to normalize a vector
void normalize(Vertex* v) {
    double length = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
    v->x /= length;
    v->y /= length;
    v->z /= length;
}

double dot_product(Vertex* v1, Vertex* v2) {
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

