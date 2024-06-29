// compile using    gcc -o raytrace main.c misc/vars.c misc/func.c -lm

#include "misc/vars.h"
#include "misc/func.h"




int main() {
    double t0 = time(NULL);

    // load vertices
    process_indicator("loading vertices", "####");
    Vertex *vertices;
    unsigned int num_thetas = 0;
    unsigned int num_phis = 0;
    read_vertices_2d(&vertices, "./vertices.txt", &num_thetas, &num_phis);
    double dphi = (double)2*PI/num_phis;


    // assign surfaces - THIS ONLY WORKS FOR SEMI-SPHERES: need to change num_surfaces otherwise
    process_indicator("assigning surfaces", "####");
    Vertex centre = {0, 0, 0};
    Surface *surfaces;
    unsigned int num_surfaces = (num_thetas - 2) * 2 * num_phis + num_phis;
    printf("Total Number of Surfaces: %i \n", num_surfaces);
    assign_surfaces(&surfaces, vertices, num_thetas, num_phis);

    // move and mirror droplet
    process_indicator("moving and mirroring droplet", "####");
    centre.z = d1;
    transform_droplet(&vertices, d1, num_thetas, num_phis);
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



    unsigned int num_beams = lens_grid_num_phi*lens_grid_num_r*circle_grid_num_r*circle_grid_num_phi;
    printf("Total Number of Beams: %i \n", num_beams);
    Beam *beams;
    beams = malloc(num_beams * sizeof(Beam));
    process_indicator("malloced beam array", "####");
    int rrc,ppc,rrl,ppl;
    unsigned int ibeam = 0;
    for (rrc = 0; rrc < circle_grid_num_r; rrc++) {
        for (ppc = 0; ppc < circle_grid_num_phi; ppc++) {
            for (rrl = 0; rrl < lens_grid_num_r; rrl++) {
                for (ppl = 0; ppl < lens_grid_num_phi; ppl++) {
                    double rc = circle_radius/(circle_grid_num_r+1)*((double) rrc+0.5);
                    double phic = (double)ppc/circle_grid_num_phi*2*PI;
                    // start the beams at different places at z = d1+d2
                    Vertex circle = {rc*cos(phic),rc*sin(phic),d1+d2};
                    beams[ibeam].length = 0;
                    beams[ibeam].x = circle.x;
                    beams[ibeam].y = circle.y;
                    beams[ibeam].z = circle.z;
                    
                    double rl = lens_radius/(lens_grid_num_r+1)*rrl;
                    double phil = (double)ppl/lens_grid_num_phi*2*PI;
                    // aim at many different points on the lens
                    Vertex lens = {rl*cos(phil),rl*sin(phil),d1};

                    Vertex diff = {lens.x-circle.x, lens.y-circle.y, lens.z-circle.z};

                    normalize(&diff);

                    beams[ibeam].vx = diff.x;
                    beams[ibeam].vy = diff.y;
                    beams[ibeam].vz = diff.z;

                    ibeam++;
                }
            }
        }
    }        

    process_indicator("initialized beams", "####");

    // print_beams(beams, num_beams);

    Vertex normal1 = {0., 0., -1.};

    double counter = 0; 
    char arrow[ARROW_LENGTH + 1]; // +1 for the null terminator   

    // MULTITHREADING IS GOOD AT SOMETHING AROUND A FEW MILLION RAYS ON MY PC
    // process_indicator("multithreading...", "####");
    // int num_procs = omp_get_num_procs();
    // omp_set_num_threads(num_procs);

    // #pragma omp parallel for schedule(static) 
    for(ibeam = 0; ibeam < num_beams; ibeam++) {
        // define variables
        Vertex intersection;
        int intersection_surface_index;

        if(ibeam == 0) process_indicator("multithreading set up!", "####");
        fprintf(outfile, "%0.20f %0.20f %0.20f , ", beams[ibeam].x, beams[ibeam].y, beams[ibeam].z);

        // first see where we intersect with the lens
        intersection_surface_index = 0;

        // one does not need to look at all surfaces!!! the surfaces are bound by the lens radius and the angle phi!
        double current_phi = (double)((ibeam % circle_grid_num_r) / (lens_grid_num_r * lens_grid_num_phi)) / circle_grid_num_phi * 2 * PI;

        // the current angle falls in the angle category 
        int current_phi_lens = (int)current_phi/dphi;

        int index1 = 0;
        int index2 = 0;
        for(int ss = 0; ss < loop_upper_bound; ss++) {
            index1 = (int)(ss==0) * 0 + (int)(ss!=0) * (int)(2*((double)ss-0.5));
            index2 = (int)(ss==0) * 0 + (int)(ss!=0) * (int)(2*((double)ss-0.5)) + 1;
            intersection_surface_index += index1*(int)intersect_triangle(&beams[ibeam], &surfaces[index1], &intersection);
            intersection_surface_index += index2*(int)intersect_triangle(&beams[ibeam], &surfaces[index2], &intersection);
        }
        
        // technically no if statement needed because no internal reflections considered
        Vertex beam_before = {beams[ibeam].x, beams[ibeam].y, beams[ibeam].z};
        beams[ibeam].x = intersection.x;
        beams[ibeam].y = intersection.y;
        beams[ibeam].z = intersection.z;
        Vertex beam_after = {beams[ibeam].x, beams[ibeam].y, beams[ibeam].z};
        Vertex beam_difference = {beam_after.x-beam_before.x, beam_after.y-beam_before.y, beam_after.z-beam_before.z};
        beams[ibeam].length += sqrt(dot_product(&beam_difference, &beam_difference));
        

        fprintf(outfile, "%0.20f %0.20f %0.20f , ", beams[ibeam].x, beams[ibeam].y, beams[ibeam].z);

        //printf("%0.20f %0.20f %0.20f \n", normals[intersection_surface_index].x, normals[intersection_surface_index].y, normals[intersection_surface_index].z);
        //printf("%0.20f %0.20f %0.20f \n", beams[ibeam].vx, beams[ibeam].vy, beams[ibeam].vz);

        beams[ibeam] = refract(beams[ibeam], normals[intersection_surface_index], nAir, nWater);


        move_beam(&beams[ibeam], (beams[ibeam].z-d1)/fabs(beams[ibeam].vz));
        beams[ibeam].length += (beams[ibeam].z-d1)/fabs(beams[ibeam].vz);

        fprintf(outfile, "%0.20f %0.20f %0.20f , ", beams[ibeam].x, beams[ibeam].y, beams[ibeam].z);

        beams[ibeam] = refract(beams[ibeam], normal1, nWater, nAir);

        move_beam(&beams[ibeam], d1/(fabs(beams[ibeam].vz)));
        beams[ibeam].length += d1/(fabs(beams[ibeam].vz));

        fprintf(outfile, "%0.20f %0.20f %0.20f , ", beams[ibeam].x, beams[ibeam].y, beams[ibeam].z);

        fprintf(outfile, "\n");


        // this is ugly though it works
        counter+=(double)ARROW_LENGTH/num_beams;
        for(int ii = 0; ii < (int)counter; ii++) {
            arrow[ii] = '-';
        }
        
        arrow[(int)counter+1] = '\0'; // Set the next character to null terminator
        printf("\r[%s>]", arrow); // Print the arrow with a carriage return
        fflush(stdout); // Flush the output buffer
    }

    printf("\n");

    process_indicator("finished raytracing", "####");

    double complex evalgrid[circle_grid_num_r];
    for(int ii = 0; ii < circle_grid_num_r; ii++) {
        evalgrid[ii] = 0.0 + 0.0*I;
    }

    for(int ii = 0; ii < num_beams; ii++) {
        int r = ii % circle_grid_num_r;
        // calculate phase based on distance travelled
        double phase = 2.0 * M_PI * beams[ii].length / wavelength;
        // calculate intensity based on where the beam ended up
        double intensity = 0.0f;
        if(beams[ii].x*beams[ii].x + beams[ii].y*beams[ii].y < dot_radius*dot_radius && beams[ii].z < d1) intensity = 1.0f;
        // add the contribution of this beam to the grid point
        evalgrid[r] += intensity * cexp(I * phase);
    }

    process_indicator("evaluated beams", "####");

    double intensity_grid[circle_grid_num_r];

    for(int ii = 0; ii < circle_grid_num_r; ii++) {
        intensity_grid[ii] = pow(cabs(evalgrid[ii]), 2);
    }

    FILE *file = fopen("intensity_grid.txt", "w");
    if (file == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for(int ii = 0; ii < circle_grid_num_r; ii++) {
        fprintf(file, "%0.20f %f \n", (double)ii/circle_grid_num_r*circle_radius, intensity_grid[ii]);
    }

    fclose(file);

    process_indicator("wrote intensities to file", "####");

    

    

    fclose(outfile);

    free(beams);
    free(vertices);
    free(surfaces);

    printf("total runtime: %i s \n ", (int)(time(NULL)-t0));
}