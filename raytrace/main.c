// compile using    gcc -o raytrace main.c misc/vars.c misc/func.c -lm

#include "misc/vars.h"
#include "misc/func.h"




int main() {
    theta_max = atan(0.0010382763/d1);
    // make phi resolution equivalent to theta resolution
    phi_res = (phi_max - phi_min) / (theta_max - theta_min) * theta_res;
    int num_total_beams = theta_res*phi_res;
    printf("Number of total beams: \%i \n", num_total_beams);

    // load vertices
    process_indicator("loading vertices", "####");
    Vertex *vertices;
    unsigned int num_thetas = 0;
    unsigned int num_phis = 0;
    read_vertices(&vertices, "./vertices.txt", &num_thetas, &num_phis);


    // assign surfaces - THIS ONLY WORKS FOR SEMI-SPHERES: need to change num_surfaces otherwise
    process_indicator("assigning surfaces", "####");
    Vertex centre = {0, 0, 0};
    Surface *surfaces;
    unsigned int num_surfaces = (num_thetas - 2) * 2 * num_phis + num_phis;
    assign_surfaces(&surfaces, vertices, num_thetas, num_phis);

    // move and mirror droplet
    process_indicator("moving and mirroring droplet", "####");
    centre.z = d1;
    transform_mirror_droplet(&vertices, d1, num_thetas, num_phis);
    print_surfaces_to_file(surfaces, num_surfaces, "./surfaces.txt", "w");

    // calculate and save surface normals
    process_indicator("calculating surface normals", "####");
    Vertex *normals = malloc(num_surfaces * sizeof(Vertex));
    calculate_surface_normals(surfaces, normals, num_surfaces, &centre);
    export_surface_normals(surfaces, normals, num_surfaces);

    // empty evaluation.txt
    process_indicator("clearing evaluation.txt", "####");
    FILE *fptr = fopen("evaluation.txt", "w");
    fclose(fptr);

    FILE *outfile = fopen("rays.txt", "w");

    process_indicator("entering main loop", "####");
    // raytrace beam over all angles
    // this code does not assume radial symmetry! I should later implement one that does
    Beam beam;
    Vertex intersection;
    int intersection_surface_index;
    bool intersect;
    for(int tt = 0; tt < theta_res; tt++) {   
        for(int pp = 0; pp < phi_res; pp++) {
            if(theta_min == 0 && tt == 0 && pp > 0)   break;      // no need to calculate the point with theta = 0 more than once
            
            initialize_beam(&beam, tt, pp);
            fprintf(fptr, "%f %f %f, ", beam.x, beam.y, beam.z);

            intersect = true;
            while(intersect == true) {
                // check for intersections
                for(int ii = 0; ii < num_surfaces; ii++) {
                    intersect = intersect_triangle(&beam, &surfaces[ii], &intersection);
                    if(intersect == true) {
                        intersection_surface_index = ii;
                        // printf("INTERSECT at %0.10f %0.10f %0.10f \n", intersection.x, intersection.y, intersection.z);
                        break;
                    }
                }

                // if intersect is true
                if(intersect == true) {
                    beam.x = intersection.x;
                    beam.y = intersection.y;
                    beam.z = intersection.z;

                    refract(&beam, &normals[intersection_surface_index]);

                } else {
                    if(beam.vz < 0) {
                        printf("ERROR! BEAM IN NEGATIVE DIRECTION! \n");

                    } else {
                        move_beam(&beam, (d1-beam.z)/beam.vz);

                    }

                }
                fprintf(fptr, "%f %f %f, ", beam.x, beam.y, beam.z);
                    
            }

            move_beam(&beam, d2/beam.vz);

            fprintf(fptr, "%f %f %f, ", beam.x, beam.y, beam.z);

            fprintf(fptr, "\n");

            // if necessary, refract

            // move until d1 (+ d2)

            // // glass refraction

            // // move by dd

            // // glass refraction

            // // move until d1 + d2 

            // move_beam(&beam, d2/beam.vz);

            
            
            Colour colour = evaluate_beam_circle(&beam, tan(theta_max/1.2)*(d1+d2));


            print_evaluation(theta_min + tt * (theta_max - theta_min) / theta_res, phi_min + pp * (phi_max - phi_min) / phi_res, colour, "evaluation.txt");
        }
    }

    fclose(outfile);

    free(vertices);
    free(surfaces);
}