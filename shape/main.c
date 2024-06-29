#include "misc/vars.h"
#include "misc/func.h"




int main() {
    // initiate 

    srand(time(NULL));
    const unsigned int n = 15;  // the lowest this code will run at is n = 2 (i think)
    Vertex* vertices;
    //Droplet droplet = create_spherical_droplet(2*n, n+1, 0.5f, &vertices);
    Droplet droplet = create_half_droplet(2*n, n/2+1, 0.001f, &vertices);
    printf("sigmas: %f %f \n", droplet.sigma_to_air, droplet.sigma_to_surface);
    normalize_volume(&droplet);
    //print_surfaces_to_file(&droplet, "surfaces.txt");
    Vertex com = calculate_center_of_mass(&droplet);
    printf("Center of mass: x=%f, y=%f, z=%f \n", com.x, com.y, com.z);
    printf("Potential energy: %f \n", PE(&droplet));
    printf("Surface energy: %f \n", SE(&droplet));

    // calculate


    unsigned int mc_timesteps = 40000;
    unsigned int write_period = 2000; 

    print_surfaces_to_file(&droplet, "surfaces.txt", "w");
    for(int ii = 0; ii < mc_timesteps; ii++) {
        mc_timestep_hemi(&droplet);
        if(ii%write_period == 0) {
            print_surfaces_to_file(&droplet, "surfaces.txt", "a");
        }
    }    

    print_specs_to_file(&droplet, "../raytrace/vertices.txt", "w");
    print_vertices_to_file(&droplet, "../raytrace/vertices.txt", "a");
    

    // destroy

    process_indicator("freeing vertices", "####");
    free(vertices);
    process_indicator("destroy droplet", "####");
    destroy_droplet(&droplet);
    return 0;
}