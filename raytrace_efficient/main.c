// compile using    gcc -o raytrace main.c misc/vars.c misc/func.c -lm

#include "misc/vars.h"
#include "misc/func.h"




int main() {
    int num_total_beams = theta_res*phi_res;
    printf("Number of total beams: %i \n", num_total_beams);

    // load vertices
    process_indicator("loading vertices", "####");
    Vertex *vertices;
    unsigned int num_thetas = 0;
    unsigned int num_phis = 0;
    read_vertices_2d(&vertices, "./vertices.txt", &num_thetas, &num_phis);
    double dphi = (double)2*PI/num_phis;
    printf("Number of thetas and phis: %i %i \n", num_thetas, num_phis);


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

    printf("max radius: %0.20f \n", lens_radius);

    // we want  to find the maximum r we have to look at surfaces for!
    int max_theta = 0;
    for(max_theta = 0; max_theta < num_thetas; max_theta ++) {
        int index = (int)(max_theta==0) * 0 + (int)(max_theta!=0) * (int)(2*((double)max_theta-0.5));
        if(surfaces[index].v3->x > lens_radius) break;
    }
    max_theta++;

    int loop_upper_bound = 2 * max_theta;


    FILE *outfile = fopen("rays.txt", "w");

    process_indicator("entering main loop", "####");
    // raytrace beam over all angles
    Beam beam;
    Vertex intersection;
    int intersection_surface_index = 0;
    unsigned int ibeam = 0;
    for(int tt = 0; tt < theta_res; tt++) {
        //for(int pp = 0; pp < phi_res; pp++) {
        for(int pp = 0; pp < 1; pp++) {
            if(theta_min == 0 && tt == 0 && pp > 0)   break;      // no need to calculate the point with theta = 0 more than once
            initialize_beam(&beam, tt, pp);
            fprintf(fptr, "%f %f %f, ", beam.x, beam.y, beam.z);

            // check for intersections
            // one does not need to look at all surfaces!!! the surfaces are bound by the lens radius and the angle phi!
            double current_phi = (double)(phi_min + (phi_max-phi_min)/phi_res * pp);

            // the current angle falls in the angle category 
            int current_phi_lens = (int)current_phi/dphi;


            int index1 = 0;
            int index2 = 0;
            intersection_surface_index = 0;
            for(int ss = 0; ss < loop_upper_bound; ss++) {
                index1 = (int)(ss==0) * 0 + (int)(ss!=0) * (int)(2*((double)ss-0.5));
                index2 = (int)(ss==0) * 0 + (int)(ss!=0) * (int)(2*((double)ss-0.5)) + 1;
                intersection_surface_index += index1*(int)intersect_triangle(&beam, &surfaces[index1], &intersection);
                intersection_surface_index += index2*(int)intersect_triangle(&beam, &surfaces[index2], &intersection);
            }


        
            beam.x = intersection.x;
            beam.y = intersection.y;
            beam.z = intersection.z;


            refract(&beam, &normals[intersection_surface_index]);

            
            if(beam.vz < 0) {
                printf("ERROR! BEAM IN NEGATIVE DIRECTION! \n");

            } else {
                move_beam(&beam, (d1-beam.z)/beam.vz);
            }

            fprintf(fptr, "%f %f %f, ", beam.x, beam.y, beam.z);
            
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

            
            
            Colour colour = evaluate_beam_lines(&beam);


            print_evaluation(theta_min + tt * (theta_max - theta_min) / theta_res, phi_min + pp * (phi_max - phi_min) / phi_res, colour, "evaluation.txt");
        
            ibeam ++;
        }
    }

    fclose(outfile);

    free(vertices);
    free(surfaces);
}