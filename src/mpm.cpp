#include "simulation.hpp"
#include <getopt.h>

int main(int argc, char* argv[]) {
     Simulation sim;
    
    // parsing management
    option longopts[] = {
        {"config", optional_argument, NULL, 'c'}, 
        {"height", required_argument, NULL, 'H'},
        {"help", optional_argument, NULL, 'h'}, {0}};

     
      while (1) {
        const int opt = getopt_long(argc, argv, "chH:", longopts, 0);

        if (opt == -1) {
            break;
        }

        switch (opt) {
            case 'c': // config file
                sim.config_loading = true;
                sim.config_file = optarg;
            case 'H': // heigth file
                // if(argc==2) {
                //     std::cout << "Please indicate the heigth file name" << std::endl;
                //     abort(); 
                // }
                std::cout << optarg << std::endl;
                sim.height_loading = true;
                sim.height_file = optarg;
                break;
            case 'h': // display help
                std::cout << "usuage is \n --config (-c) : for enabling the reading of config file" << std::endl;
                std::cout << "--height (-H) : for enabling the reading of a height file (format in readme.txt)" << std::endl;
                abort();
            default:
                std::cout << "option not defined" << std::endl;
                abort();
        } 
    }

      // sim.sim_name = "MCC_p0_3_E1_M3_dx_10_PIC90";
    //   sim.sim_name = "MCC_safeFa05_NApic_real_15000_17_01";
    //   sim.sim_name = "Elastic_safeFa05_dx025_LARGE_centered";
        // sim.sim_name = "Napic_flip09_safeFa05";
      sim.sim_name = "PST_USL_PIC09";
    //   sim.directory = "../../../Simulation/Papier/Square/";
    sim.directory = "/media/guillet/HDD1/Guillet/Simulation/PST/";
      //   //For benchmark number 1
      //   sim.sim_name = std::to_string( (int) sim.L1 ) + "Cohesion_WL=" + std::to_string( (int) sim.Initial_cohesion_WL );


      sim.sim_type = "release"; // only "flow" and "release" choices
      sim.compressibility = false; // boolean
      sim.parallelised = true;
      sim.plasticity_model = "No Plasticity"; // for now only "Drucker-Prager", "No Plasticity" and "MCC"
      sim.saving_type = "Light"; //only Light, ALL, PST, CrossSlope
      sim.poisson_disk = false;

      debug("Plasticity Model = ", sim.plasticity_model, "\n");
      debug("Simulation Type  = ", sim.sim_type, "\n");

      T material_density = 250.;
      T slope_angle = 45./180.*M_PI;
      sim.end_frame = 240;
    //   sim.end_frame = 1200.;
      sim.frame_dt = 1.0 / 240.0;
    // sim.frame_dt = 1./480.;
      // sim.dx = 1./(20.*cos(slope_angle));
      // sim.dx = cos(slope_angle)/10.;
      sim.dx = 1./50.;
    //   sim.dx = 0.25;
      
      sim.loading_time = 0.25; // must be higher than the simulation time
      //L2 is in the Y direction qnd L1 in the X direction
    //   sim.velocity_PST = 3.;
      sim.velocity_PST = 6.;
      
      // T slope_angle = 45./180.*M_PI;
      
      sim.L1 = 0.3;
      // sim.L2 = 20./cos(slope_angle);
      sim.L2 = 40.;

      sim.gravity = 9.80665;
      sim.cfl = 0.5;
    //   sim.flip_ratio = 0.;
    sim.flip_ratio = 0.9;
      sim.basal_frictional_angle = 27./180.*M_PI;
      sim.fluidization = 0.;

      sim.n_threads = 6;
      sim.rho_water = 1000.;
      sim.phi = 30./180.*M_PI; // internal frictional angle needs to be non constant
      sim.psi = 0./180.*M_PI; // dilatancy angle
      sim.cohesion = 1e8; // cohesion parameter
      sim.nu = 0.3;
      sim.E = 10.e6;
      sim.Hp = 0.;//sim.E;
      sim.safety_factor = 0.5;
      sim.void_ratio = 0.8;
      sim.mu = 0.; // dry sand = 0.

      sim.Shear_strength = 1450.;
      // A changer

      sim.Initial_cohesion_WL = sim.Shear_strength - material_density*sim.gravity*0.5*cos(slope_angle)*tan(sim.basal_frictional_angle);

      sim.Thickness_WL = 0.125;
      //sim.Shear_modulus_WL = 2.5e5;
      // sim.Shear_modulus_WL = 1e6/(2*(1 + sim.nu));
      sim.Shear_modulus_WL = 1e6/(2*(1.3));
      sim.Stiffness_WL = sim.Shear_modulus_WL/sim.Thickness_WL;

      //MCC quantities:
    //   sim.beta = 0.3;
      sim.beta = 0.2;
      //sim.M = 30./180.*M_PI;
    //   sim.M = 2.7;
      sim.M = 1.7;
      sim.p0_Initial = 8000.;
    // sim.p0_Initial = 12000;
      sim.threshold_MCC = 1e-3;

    //   unsigned int Nloop = std::round(sim.L/sim.dx);
    //   debug("Nloop           = ", Nloop);

      unsigned int NloopX = std::round(sim.L1/sim.dx);
      unsigned int NloopY = std::round(sim.L2/sim.dx);
      debug("NloopX           = ", NloopX, "\n");
      debug("NloopY           = ", NloopY, "\n");

      const unsigned int ppc = 4;
      sim.Nx = NloopX;
      sim.Ny = NloopY;
      sim.Np = NloopX * NloopY * ppc;

      sim.initialize(material_density);

      debug("Wave speed      = ", sim.wave_speed);
      debug("dt_max          = ", sim.dt_max);
      debug("Np              = ", sim.Np);

      sim.initialize_particles(/* number of loops in x */ NloopX,
                             /* number of loops in x */ NloopY,
                             /* nb of part. per cell*/ ppc,
                             /* material density    */ material_density,
                             /* slope angle in y */ slope_angle);


    
    sim.simulate();

	return 0;
}
