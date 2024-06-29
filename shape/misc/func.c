#include "func.h"
#include "vars.h"




// ------------------------------------------- PHYSICAL FUNCTIONS


double calc_volume(Droplet *droplet) {
    double total_volume = 0.0;
    for (int i = 0; i < droplet->num_surfaces; i++) {
        Surface* surface = &(droplet->surfaces[i]);

        // calculate the volume of the tetrahedron
        double volume = vertex_span_volume(*surface->v1, *surface->v2, *surface->v3);

        // add the volume of the tetrahedron to the total volume
        total_volume += volume;
    }
    return total_volume;
}


Vertex calculate_center_of_mass(Droplet* droplet) {
    Vertex center = {0.0, 0.0, 0.0};
    float total_volume = 0.0;

    // iterate over all surfaces in the droplet
    for (int i = 0; i < droplet->num_surfaces; i++) {
        Surface* surface = &(droplet->surfaces[i]);

        // calculate the center of mass of the tetrahedron
        Vertex tetrahedron_center = vertex_scalar_mult(vertex_add(vertex_add(*surface->v1,*surface->v2),*surface->v1),1.0/4.0);

        // calculate the volume of the tetrahedron
        double volume = vertex_span_volume(*surface->v1, *surface->v2, *surface->v3);

        // add the center of mass of the tetrahedron, weighted by its volume, to the total
        center = vertex_add(center, vertex_scalar_mult(tetrahedron_center, volume));

        // add the volume of the tetrahedron to the total volume
        total_volume += volume;
    }
    // printf("total volume: %f \n", total_volume);

    // divide by the total volume to get the weighted average position
    center = vertex_scalar_mult(center, 1.0/total_volume);

    return center;
}


double PE (Droplet *droplet) {
    return droplet->density*calc_volume(droplet)*g*calculate_center_of_mass(droplet).z;
}


double SE(Droplet* droplet) {
    double total_surface_energy = 0.0;
    double total_area = 0.0;

    // iterate over all surfaces in the droplet
    for (int i = 0; i < droplet->num_surfaces; i++) {
        Surface* surface = &(droplet->surfaces[i]);

        // calculate the vectors representing two sides of the triangle
        Vertex AB = vertex_subtract(*surface->v2, *surface->v1);
        Vertex AC = vertex_subtract(*surface->v3, *surface->v1);

        // calculate the magnitude of the cross product
        double area = vertex_span_area(AB, AC);

        // add to total area
        total_area += area;
    }
    total_surface_energy = droplet->sigma_to_air * total_area;

    return total_surface_energy;
}


double SE_bottom(Droplet* droplet) {
    double total_surface_energy = 0.0;
    double total_area = 0.0;

    // iterate over the bottom surfaces of the droplet
    for(int ii = 0; ii < droplet->num_phis; ii++) {
        Vertex v1 = droplet->ptr_vertices[0][droplet->num_phis - ii];
        Vertex v2 = droplet->ptr_vertices[0][droplet->num_phis - (ii + 1) % droplet->num_phis];

        total_area += vertex_span_area(v1, v2);
    }

    total_surface_energy = total_area * droplet->sigma_to_surface;

    return total_surface_energy;
}

// IMPORTANT!!! RADIAL SYMMETRY WAS ASSUMED EXPLICITLY HERE
// This function should be built for efficiency - although this will probably not work well here

// This function DOES NOT RESCALE VOLUME
double delta_PE_hemi(Droplet* droplet, unsigned int theta, double scale_factor) {
    // sadly, this will ahve to scale with n^2, as the volume needs to be normalized...
    // i could make it more efficient (O(n)) by making the integrals one-dimensional!! and only looking at z!
    double PE1 = PE(droplet);

    scale_row_hemi(droplet, theta, scale_factor);

    double PE2 = PE(droplet);

    scale_row_hemi(droplet, theta, 1.0/scale_factor);

    return PE2 - PE1;
}


// IMPORTANT!!! RADIAL SYMMETRY WAS ASSUMED EXPLICITLY HERE
// This function should be built for efficiency

// this function DOES NOT RESCALE VOLUME
double delta_SE_hemi(Droplet* droplet, unsigned int theta, double scale_factor) {
    double delta_SE = 0;

    if(theta == 0)  {
        // calculate old surface energy

        Vertex v0 = droplet->ptr_vertices[0][0];
        Vertex v1 = droplet->ptr_vertices[0][1];
        Vertex v2 = droplet->ptr_vertices[0][2];
        
        Vertex AB = vertex_subtract(v1,v0);
        Vertex AC = vertex_subtract(v2,v0);

        double area = vertex_span_area(AB,AC) * droplet->num_phis;
        double SE1 = area * droplet->sigma_to_air;


        // move top vertex 

        v0 = vertex_scalar_mult(v0, scale_factor);


        // calculate new surface energy

        AB = vertex_subtract(v1,v0);
        AC = vertex_subtract(v2,v0);

        area = vertex_span_area(AB,AC) * droplet->num_phis;
        double SE2 = area * droplet->sigma_to_air;

        delta_SE = SE2 - SE1;

    } else if (theta == droplet->num_thetas-1) {
        // calculate old surface energy

        Vertex v0 = droplet->ptr_vertices[0][droplet->num_vertices - 1];
        Vertex v1 = droplet->ptr_vertices[0][droplet->num_vertices - 2];
        Vertex v2 = droplet->ptr_vertices[0][droplet->num_vertices - 1 - droplet->num_phis];
        Vertex v3 = droplet->ptr_vertices[0][droplet->num_vertices - 2 - droplet->num_phis];

        Vertex AB = vertex_subtract(v1,v0);
        Vertex AC = vertex_subtract(v2,v0);

        Vertex DB = vertex_subtract(v1,v3);
        Vertex DC = vertex_subtract(v2,v3);

        double area = (vertex_span_area(AB,AC) + vertex_span_area(DB,DC)) * droplet->num_phis;
        double SE1 = area * droplet->sigma_to_air;

        // the bottom surface

        double area_bot = vertex_span_area(v0,v1) * droplet->num_phis;
        SE1 += area_bot * droplet->sigma_to_surface;


        // move vertices

        v0 = vertex_scalar_mult(v0, scale_factor);
        v1 = vertex_scalar_mult(v1, scale_factor);


        // calculate new surface energy

        AB = vertex_subtract(v1,v0);
        AC = vertex_subtract(v2,v0);

        DB = vertex_subtract(v1,v3);
        DC = vertex_subtract(v2,v3);

        
        area = (vertex_span_area(AB,AC) + vertex_span_area(DB,DC)) * droplet->num_phis;

        double SE2 = area * droplet->sigma_to_air;


        // the bottom surface

        area_bot = vertex_span_area(v0,v1) * droplet->num_phis;
        SE2 += area_bot * droplet->sigma_to_surface;


        // difference

        delta_SE = SE2 - SE1;

    } else if (theta == 1){
        // old surface energy

        unsigned int offset = 1;

        Vertex v0 = droplet->ptr_vertices[0][offset];
        Vertex v1 = droplet->ptr_vertices[0][offset + 1];
        Vertex v2 = droplet->ptr_vertices[0][offset + droplet->num_phis];
        Vertex v3 = droplet->ptr_vertices[0][offset + 1 + droplet->num_phis];
        Vertex v4 = droplet->ptr_vertices[0][0];


        Vertex AB = vertex_subtract(v1,v0);
        Vertex AD = vertex_subtract(v3,v0);
        Vertex AC = vertex_subtract(v2,v0);
        Vertex AE = vertex_subtract(v4,v0);

        double area = (vertex_span_area(AB,AD) + vertex_span_area(AC,AD) + vertex_span_area(AB,AE)) * droplet->num_phis;
        double SE1 = area * droplet->sigma_to_air;


        // move vertices

        v0 = vertex_scalar_mult(v0, scale_factor);
        v1 = vertex_scalar_mult(v1, scale_factor);


        // calculate new surface energy

        AB = vertex_subtract(v1,v0);
        AD = vertex_subtract(v3,v0);
        AC = vertex_subtract(v2,v0);
        AE = vertex_subtract(v4,v0);

        area = (vertex_span_area(AB,AD) + vertex_span_area(AC,AD) + vertex_span_area(AB,AE)) * droplet->num_phis;
        double SE2 = area * droplet->sigma_to_air;


        // difference

        delta_SE = SE2 - SE1;
        
    } else if (theta < 0 || theta >= droplet->num_thetas) {
        printf("ERROR: trying to scale row that is out of bounds: theta = %i, num_thetas = %i", theta, droplet->num_thetas);
    } else {
        // could be simplified via trapezoids??

        // old surface energy

        unsigned int offset = 1 + droplet->num_phis * (theta - 1);

        Vertex v0 = droplet->ptr_vertices[0][offset];
        Vertex v1 = droplet->ptr_vertices[0][offset + 1];
        Vertex v2 = droplet->ptr_vertices[0][offset + droplet->num_phis];
        Vertex v3 = droplet->ptr_vertices[0][offset + 1 + droplet->num_phis];
        Vertex v4 = droplet->ptr_vertices[0][offset - droplet->num_phis];
        Vertex v5 = droplet->ptr_vertices[0][offset + 1 - droplet->num_phis];


        Vertex AB = vertex_subtract(v1,v0);
        Vertex AD = vertex_subtract(v3,v0);
        Vertex AC = vertex_subtract(v2,v0);
        Vertex AE = vertex_subtract(v4,v0);
        Vertex DF = vertex_subtract(v5,v3);
        Vertex DB = vertex_subtract(v1,v3);

        double area = (vertex_span_area(AB,AD) + vertex_span_area(AC,AD) + vertex_span_area(AB,AE) + vertex_span_area(DF,DB)) * droplet->num_phis;
        double SE1 = area * droplet->sigma_to_air;


        // move vertices

        v0 = vertex_scalar_mult(v0, scale_factor);
        v1 = vertex_scalar_mult(v1, scale_factor);


        // calculate new surface energy

        AB = vertex_subtract(v1,v0);
        AD = vertex_subtract(v3,v0);
        AC = vertex_subtract(v2,v0);
        AE = vertex_subtract(v4,v0);
        DF = vertex_subtract(v5,v3);
        DB = vertex_subtract(v1,v3);

        area = (vertex_span_area(AB,AD) + vertex_span_area(AC,AD) + vertex_span_area(AB,AE) + vertex_span_area(DF,DB)) * droplet->num_phis;
        double SE2 = area * droplet->sigma_to_air;
        

        // difference

        delta_SE = SE2 - SE1;
    }
    return delta_SE;
}














// ------------------------------------------- (RADIALLY SYMMETRIC) MONTE CARLO FUNCTIONS
// - take a random row of constant theta and displace radially by a certain factor. 
// - normalize the volume to V0

void scale_droplet(Droplet *droplet, double scale_factor) {
    for(int ii = 0; ii < droplet->num_vertices; ii++) {
        droplet->ptr_vertices[0][ii] = vertex_scalar_mult(droplet->ptr_vertices[0][ii], scale_factor);
    }
}

void normalize_volume(Droplet *droplet) {
    double V = calc_volume(droplet);
    double V_ratio = droplet->volume0/V;
    double scale_factor = pow(V_ratio, 1.0f/3.0f);
    scale_droplet(droplet, scale_factor);
}

void scale_row_spherical(Droplet *droplet, unsigned int theta, double scale_factor) {
    if(theta == 0)  {
        droplet->ptr_vertices[0][0] = vertex_scalar_mult(droplet->ptr_vertices[0][0], scale_factor);
    } else if (theta == droplet->num_thetas-1) {
        droplet->ptr_vertices[0][droplet->num_vertices-1] = vertex_scalar_mult(droplet->ptr_vertices[0][droplet->num_vertices-1], scale_factor);
    } else if (theta < 0 || theta >= droplet->num_thetas) {
        printf("ERROR: trying to scale row that is out of bounds: theta = %i, num_thetas = %i", theta, droplet->num_thetas);
    } else {
        unsigned int offset = 1 + droplet->num_phis * (theta - 1);
        for(int ii = 0; ii < droplet->num_phis; ii++) {
            droplet->ptr_vertices[0][offset + ii] = vertex_scalar_mult(droplet->ptr_vertices[0][offset + ii], scale_factor);
        }
    }
    // don't forget to normalize volume after
    normalize_volume(droplet);
}

double mc_etot_spherical(Droplet *droplet) {
    return SE(droplet) + PE(droplet);
}

void try_scale_random_row_spherical(Droplet *droplet) {
    double e0 = mc_etot_spherical(droplet);

    unsigned int random_row = rand() % droplet->num_thetas;
    double random_scale_factor = 1.0f + epsilon * (double)rand() / RAND_MAX;
    if((double)rand() / RAND_MAX < 0.5f) {
        random_scale_factor = 1/random_scale_factor;
    }
    scale_row_spherical(droplet, random_row, random_scale_factor);

    double e1 = mc_etot_spherical(droplet);

    // "zero temperature" mc acceptance rates
    if(e1 > e0) {
        scale_row_spherical(droplet, random_row, 1.0/random_scale_factor);
        //printf("scale denied \n");
    } else {
        //printf("scale accepted \n");
    }
}

void mc_timestep_spherical(Droplet *droplet) {
    for(int ii = 0; ii < droplet->num_thetas; ii++) {
        try_scale_random_row_spherical(droplet);
    }
}


void scale_row_hemi(Droplet *droplet, unsigned int theta, double scale_factor) {
    if(theta == 0)  {
        droplet->ptr_vertices[0][0] = vertex_scalar_mult(droplet->ptr_vertices[0][0], scale_factor);
    } else if (theta < 0 || theta >= droplet->num_thetas) {
        printf("ERROR: trying to scale row that is out of bounds: theta = %i, num_thetas = %i", theta, droplet->num_thetas);
    } else {
        unsigned int offset = 1 + droplet->num_phis * (theta - 1);
        for(int ii = 0; ii < droplet->num_phis; ii++) {
            droplet->ptr_vertices[0][offset + ii] = vertex_scalar_mult(droplet->ptr_vertices[0][offset + ii], scale_factor);
        }
    }
    // don't forget to normalize volume after
    normalize_volume(droplet);
}

double mc_etot_hemi(Droplet *droplet) {
    return SE(droplet) + SE_bottom(droplet);// + PE(droplet);
}

void try_scale_random_row_hemi(Droplet *droplet) {
    double e0 = mc_etot_hemi(droplet);

    unsigned int random_row = rand() % droplet->num_thetas;
    double random_scale_factor = 1.0f + epsilon * (double)rand() / RAND_MAX;
    if((double)rand() / RAND_MAX < 0.5f) {
        random_scale_factor = 1/random_scale_factor;
    }
    scale_row_hemi(droplet, random_row, random_scale_factor);

    double e1 = mc_etot_hemi(droplet);


    // "zero temperature" mc acceptance rates
    if(e1 > e0) {
        scale_row_hemi(droplet, random_row, 1.0/random_scale_factor);
        //printf("scale denied \n");
    } else {
        //printf("scale accepted \n");
    }
}

void try_scale_random_row_hemi_2(Droplet *droplet) {
    unsigned int random_row = rand() % droplet->num_thetas;
    double random_scale_factor = 1.0f + epsilon * (double)rand() / RAND_MAX;
    if((double)rand() / RAND_MAX < 0.5f) {
        random_scale_factor = 1/random_scale_factor;
    }
    double deltaSE = delta_SE_hemi(droplet, random_row, random_scale_factor);   // does not rescale volume
    double deltaPE = delta_PE_hemi(droplet, random_row, random_scale_factor);   // does not rescale volume
    double deltaV = abs(calc_volume(droplet) - droplet->volume0);
    scale_row_hemi(droplet, random_row, random_scale_factor);
    double deltaVp = abs(calc_volume(droplet) - droplet->volume0);
    scale_row_hemi(droplet, random_row, 1.0/random_scale_factor);
    double delta_E = deltaSE + deltaPE + lambda * (deltaVp - deltaV);
    printf("delta SE = %0.20f \t delta PE = %0.20f \n", deltaSE, deltaPE);
    // "zero temperature" mc acceptance rates
    if(T == 0.) {
        if(delta_E > 0.) {
            //printf("scale denied \n");
        } else {
            scale_row_hemi(droplet, random_row, random_scale_factor);
            //printf("scale accepted \n");
        }
    } else {
        if(delta_E < 0.) {
            //printf("scale denied \n");
            scale_row_hemi(droplet, random_row, random_scale_factor);
        } else {
            // maybe approximate rates?
            double rate = exp(-delta_E/T);
            if((double)rand() / RAND_MAX < rate) scale_row_hemi(droplet, random_row, random_scale_factor);
            //printf("scale accepted \n");
        }
    }
    
}

void mc_timestep_hemi(Droplet *droplet) {
    for(int ii = 0; ii < droplet->num_thetas; ii++) {
        try_scale_random_row_hemi_2(droplet);
    }
}












// ------------------------------------------- DROPLET FUNCTIONS


// for phi and theta to have about the same resolution, one should choose num_phis = 2n, num_thetas = n+1
// one can access the internal structure of the droplet via droplet.surfaces[i].v1->x
// Droplet create_spherical_droplet(int num_phis, int num_thetas, double radius, Vertex **vertices) {
//     Droplet droplet;
//     droplet.density = 1e3;
//     droplet.volume0 = 4.0f / 3.0f * PI * radius * radius * radius;
//     droplet.num_surfaces = (num_thetas - 3) * 2 * num_phis + 2 * num_phis;
//     droplet.surfaces = (Surface *)malloc(droplet.num_surfaces * sizeof(Surface));
//     droplet.ptr_vertices = vertices;
//     droplet.num_thetas = num_thetas;
//     droplet.num_phis = num_phis;
//     droplet.sigma_to_air = 0.07275;
//     droplet.sigma_to_surface = 0.04;

//     // Calculate the positions of the vertices
//     unsigned int num_vertices = (num_phis * (num_thetas - 2) + 2);
//     droplet.num_vertices = num_vertices;
//     process_indicator("mallocing vertices", "####");
//     // printf("number of vertices: %i \n", num_vertices);
//     *vertices = (Vertex *)malloc(num_vertices * sizeof(Vertex));

//     // upper pole without redundancy
//     process_indicator("assigning vertices", "####");
//     FILE *vertex_assignments = fopen("./tests/vertex_assignments.txt", "w");

//     (*vertices)[0].x = 0.0;
//     (*vertices)[0].y = 0.0;
//     (*vertices)[0].z = radius;
//     fprintf(vertex_assignments, "0 ");
//     for (int j = 0; j < num_thetas - 2; j++) {
//         for (int i = 0; i < num_phis; i++) {
//             double phi = 2 * PI * i / num_phis;
//             double theta = PI * (j+1) / (num_thetas - 1);
//             (*vertices)[j * num_phis + i + 1].x = radius * sin(theta) * cos(phi);
//             (*vertices)[j * num_phis + i + 1].y = radius * sin(theta) * sin(phi);
//             (*vertices)[j * num_phis + i + 1].z = radius * cos(theta);
//             fprintf(vertex_assignments, "%i ", j * num_phis + i + 1);
//         }
//     }
//     // lower pole without redundancy
//     (*vertices)[num_vertices - 1].x = 0.0;
//     (*vertices)[num_vertices - 1].y = 0.0;
//     (*vertices)[num_vertices - 1].z = -radius;

//     fprintf(vertex_assignments, "%i ", num_vertices - 1);


//     process_indicator("assigned vertices", "####");
//     fclose(vertex_assignments);
//     // Assign the vertices to the surfaces

//     process_indicator("assigning first surfaces", "####");
//     FILE *surface_assignments = fopen("./tests/surface_assignments.txt", "w");


//     // first the ones at the top poles
//     for (int i = 0; i < num_phis; i++) {
//         // Create one surface for each grid cell
//         droplet.surfaces[i].v1 = &(*vertices)[0];
//         droplet.surfaces[i].v2 = &(*vertices)[(i + 1)];
//         droplet.surfaces[i].v3 = &(*vertices)[(i + 1) % num_phis + 1];
//         fprintf(surface_assignments, "%i ", i);
//     }

//     process_indicator("assigning more surfaces", "####");

//     unsigned int vertex_offset = 1;
//     unsigned int surface_offset = num_phis;

//     // then the middle surfaces
//     for (int j = 0; j < num_thetas - 3; j++) {
//         for (int i = 0; i < num_phis; i++) {
//             // Create two surfaces for each grid cell
//             droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i].v1 = &(*vertices)[vertex_offset + j*num_phis + i];
//             droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i].v2 = &(*vertices)[vertex_offset + j*num_phis + ((i + 1) % num_phis)];
//             droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i].v3 = &(*vertices)[vertex_offset + (j+1)*num_phis + i];
//             fprintf(surface_assignments, "%i ", 2 * j * num_phis + surface_offset + 2 * i);

//             droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i + 1].v1 = &(*vertices)[vertex_offset + j*num_phis + ((i + 1) % num_phis)];
//             droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i + 1].v2 = &(*vertices)[vertex_offset + (j+1)*num_phis + i];
//             droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i + 1].v3 = &(*vertices)[vertex_offset + (j+1)*num_phis + ((i + 1) % num_phis)];
//             fprintf(surface_assignments, "%i ", 2 * j * num_phis + surface_offset + 2 * i + 1);
//         }
//     }

//     process_indicator("assigning even more surfaces", "####");

//     vertex_offset = (num_thetas - 3) * num_phis + 1;
//     surface_offset = (num_thetas - 3) * 2 * num_phis + num_phis;

//     for (int i = 0; i < num_phis; i++) {
//         // Create one surface for each grid cell
//         droplet.surfaces[i + surface_offset].v1 = &(*vertices)[num_vertices-1];
//         droplet.surfaces[i + surface_offset].v2 = &(*vertices)[(i + 1) + vertex_offset];
//         droplet.surfaces[i + surface_offset].v3 = &(*vertices)[(i + 1) % num_phis + 1 + vertex_offset];
//         fprintf(surface_assignments, "%i ", i + surface_offset);
//     }

//     process_indicator("assiged surfaces", "####");
//     fclose(surface_assignments);

//     return droplet;
// }


Droplet create_half_droplet(int num_phis, int num_thetas, double radius, Vertex **vertices) {
    Droplet droplet;
    droplet.density = 1e3;
    droplet.volume0 = 4.0f / 3.0f * PI * radius * radius * radius / 2;
    droplet.num_surfaces = (num_thetas - 2) * 2 * num_phis + num_phis;
    droplet.surfaces = (Surface *)malloc(droplet.num_surfaces * sizeof(Surface));
    droplet.ptr_vertices = vertices;
    droplet.num_thetas = num_thetas;
    droplet.num_phis = num_phis;
    droplet.sigma_to_air = 100;
    droplet.sigma_to_surface = 50;

    // Calculate the positions of the vertices
    unsigned int num_vertices = (num_phis * (num_thetas - 1) + 1);
    droplet.num_vertices = num_vertices;
    process_indicator("mallocing vertices", "####");
    // printf("number of vertices: %i \n", num_vertices);
    *vertices = (Vertex *)malloc(num_vertices * sizeof(Vertex));

    // upper pole without redundancy
    process_indicator("assigning vertices", "####");
    FILE *vertex_assignments = fopen("./tests/vertex_assignments.txt", "w");

    (*vertices)[0].x = 0.0;
    (*vertices)[0].y = 0.0;
    (*vertices)[0].z = radius;
    fprintf(vertex_assignments, "0 ");
    for (int j = 0; j < num_thetas - 1; j++) {
        for (int i = 0; i < num_phis; i++) {
            double phi = 2 * PI * i / num_phis;
            double theta = PI / 2.0f * (j+1) / (num_thetas - 1);
            (*vertices)[j * num_phis + i + 1].x = radius * sin(theta) * cos(phi);
            (*vertices)[j * num_phis + i + 1].y = radius * sin(theta) * sin(phi);
            (*vertices)[j * num_phis + i + 1].z = radius * cos(theta);
            fprintf(vertex_assignments, "%i ", j * num_phis + i + 1);
        }
    }


    process_indicator("assigned vertices", "####");
    fclose(vertex_assignments);
    // Assign the vertices to the surfaces

    process_indicator("assigning first surfaces", "####");
    FILE *surface_assignments = fopen("./tests/surface_assignments.txt", "w");


    // first the ones at the top poles
    for (int i = 0; i < num_phis; i++) {
        // Create one surface for each grid cell
        droplet.surfaces[i].v1 = &(*vertices)[0];
        droplet.surfaces[i].v2 = &(*vertices)[(i + 1)];
        droplet.surfaces[i].v3 = &(*vertices)[(i + 1) % num_phis + 1];
        fprintf(surface_assignments, "%i ", i);
    }

    process_indicator("assigning more surfaces", "####");

    unsigned int vertex_offset = 1;
    unsigned int surface_offset = num_phis;

    // then the middle surfaces
    for (int j = 0; j < num_thetas - 2; j++) {
        for (int i = 0; i < num_phis; i++) {
            // Create two surfaces for each grid cell
            droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i].v1 = &(*vertices)[vertex_offset + j*num_phis + i];
            droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i].v2 = &(*vertices)[vertex_offset + j*num_phis + ((i + 1) % num_phis)];
            droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i].v3 = &(*vertices)[vertex_offset + (j+1)*num_phis + i];
            fprintf(surface_assignments, "%i ", 2 * j * num_phis + surface_offset + 2 * i);

            droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i + 1].v1 = &(*vertices)[vertex_offset + j*num_phis + ((i + 1) % num_phis)];
            droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i + 1].v2 = &(*vertices)[vertex_offset + (j+1)*num_phis + i];
            droplet.surfaces[2 * j * num_phis + surface_offset + 2 * i + 1].v3 = &(*vertices)[vertex_offset + (j+1)*num_phis + ((i + 1) % num_phis)];
            fprintf(surface_assignments, "%i ", 2 * j * num_phis + surface_offset + 2 * i + 1);
        }
    }

    process_indicator("assiged surfaces", "####");
    fclose(surface_assignments);

    return droplet;
}



void destroy_droplet(Droplet *droplet) {
    // Free each surface
    process_indicator("freeing each surface","####");
    free(droplet->surfaces);
}









// ------------------------------------------- HELPER FUNCTIONS

void process_indicator(char* s, char* status){
    printf("%-*s%s\n", FIELD_WIDTH, s ,status);
}


void print_surfaces_to_file(Droplet *droplet, const char *filename, const char *type) {
    FILE *file = fopen(filename, type);
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    // printf("num surfaces: %i \n", droplet->num_surfaces);

    for (int i = 0; i < droplet->num_surfaces; i++) {
        // printf("trying to write to file: surface number %i \n", i);
        Surface *s = &droplet->surfaces[i];
        fprintf(file, "%0.10f %0.10f %0.10f, %0.10f %0.10f %0.10f, %0.10f %0.10f %0.10f\n",
                s->v1->x, s->v1->y, s->v1->z,
                s->v2->x, s->v2->y, s->v2->z,
                s->v3->x, s->v3->y, s->v3->z);
    }
    fprintf(file,"\n");

    fclose(file);
}

void print_vertices_to_file(Droplet *droplet, const char *filename, const char *type) {
    FILE *file = fopen(filename, type);
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    // printf("num surfaces: %i \n", droplet->num_surfaces);

    for (int i = 0; i < droplet->num_vertices; i++) {
        // printf("trying to write to file: surface number %i \n", i);
        Vertex *v = &droplet->ptr_vertices[0][i];
        fprintf(file, "%0.10f %0.10f %0.10f \n", v->x, v->y, v->z);
    }
    fprintf(file,"\n");

    fclose(file);
}

void print_specs_to_file(Droplet *droplet, const char *filename, const char *type) {
    FILE *file = fopen(filename, type);
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    fprintf(file, "%i %i \n", droplet->num_thetas, droplet->num_phis);
    fclose(file);
}




// ------------------------------------------- VERTEX CALC

// HAVE NOT BEEN TESTED
Vertex vertex_add(Vertex v1, Vertex v2) {
    Vertex v3 = {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
    return v3;
}

Vertex vertex_subtract(Vertex v1, Vertex v2) {
    return vertex_add(v2,vertex_scalar_mult(v1,-1.0));
}

Vertex vertex_scalar_mult(Vertex v1, double s) {
    Vertex v2 = {v1.x * s, v1.y * s, v1.z * s};
    return v2;
}

Vertex vertex_cross_product(Vertex v1, Vertex v2) {
    Vertex v3 = {v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x};
    return v3;
}

double vertex_inner_product(Vertex v1, Vertex v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double vertex_norm(Vertex v1) {
    return sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
}

double vertex_span_area(Vertex v1, Vertex v2) {
    return vertex_norm(vertex_cross_product(v1, v2))/2.0;
}

double vertex_span_volume(Vertex v1, Vertex v2, Vertex v3) {
    return fabs(vertex_inner_product(vertex_cross_product(v1,v2),v3))/6.0;         // voulume of tetrahedron is 1/6 of parrallelepiped
}

void reassign(Vertex *v1, double x, double y, double z) {
    v1->x = x;
    v1->y = y;
    v1->z = z;
}


