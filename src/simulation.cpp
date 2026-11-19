#include "simulation.hpp"

Simulation::Simulation(){
    current_time_step = 0;
    time  = 0;
    frame = 0;
    exit  = 0;

    runtime_p2g = 0;
    runtime_g2p = 0;
    runtime_g2p_1 = 0;
    runtime_g2p_2 = 0;

    // initial_k = 1.0;
    /////////// Default - will be overwritten by remesh if used ////////////
    Nx = 21; // N = (L / dx + 1)
    Ny = 21;
    grid = Grid(Nx, Ny);
    /////////////////////////////////////////////////////////////////////////
}


void Simulation::initialize(T density){

    // this function needs to be either deleted and we initialize everyting in mpm.cpp (not recommend)
    // or we rewrite it properly with all values that might change between 2 simulations

    dim = 2;

    // complete initial parameters:
    rho = density;

    wave_speed = std::sqrt(E/rho);
    // debug(" wave_speed =    ",wave_speed);
    dt_max = safety_factor * dx / wave_speed;
    // dt_max = 0.5 * dx / wave_speed;
    // dt_max = 0.01;

    // particle_area = std::pow(L, dim) / Np; // in 2D

    particles = Particles(Np);
    // std::fill( particles.A.begin(),   particles.A.end(),  std::pow(L, dim) / Np );
    std::fill( particles.A.begin(),   particles.A.end(),  (L1 * L2) / Np );

    // elastic constants
    G = E/(2.*(1.+nu));
    K = E/(2.*(1.-nu));
    G_3D = E/(2.*(1.+nu));
    K_3D = E/(3.*(1.-2.*nu));

    // Drucker Prager
    // quantities:
    eta = 6.*sin(phi)/(sqrt(3.)*(3.-sin(phi)));
    eta_g = 6.*sin(psi)/(sqrt(3.)*(3.-sin(psi)));
    // xi = 6.*cos(phi)/(sqrt(3.)*(3.+sin(phi)));
    xi = 30.;
    
    // std::cout << G << "   " << K << "   " << alpha <<"   " <<  alpha_g << "   " << xi << std::endl;

    if (sim_type == "release") {
        sim_type_int = 0;
    } else {
        sim_type_int = 1;
    }

    if (plasticity_model == "Drucker-Prager"){
        plasticity_model_int = 1;
    } else if (plasticity_model == "MCC") {
        plasticity_model_int = 2;
    } else if (plasticity_model == "No Plasticity") {
        plasticity_model_int = 0;
    } else {
        std::cout << "Invalid plasticity model" << std::endl;
    }

}
  void Simulation::initialize_particles(unsigned int const& NloopX, unsigned int const& NloopY, unsigned int const& ppc,
                                        T const& material_density, T const& slope_angle){ // slope angle to pass terrain grid object

      T K0 = nu/(1.-nu); // earth pressure

      //we srart with h initiation (or reading)
      if(height_loading) {
          read_height();
      } else {
          for (int p=0; p < Np; p++) {particles.h[p] = 0.5;}
      }

      if(poisson_disk){
          std::vector<TV2> particles_position;
          Poisson_disk fish(dx, L1, L2, ppc);
          particles_position = fish.poisson_disk_sampling();
          for(int p = 0; p < Np; p++){
            particles.x[p] = particles_position[p];
            particles.x0[p] = particles_position[p];
            }
        } else {

        std::vector<T> disp_i(ppc); disp_i[0] = 0.25; disp_i[1] = 0.75; disp_i[2] = 0.25; disp_i[3] = 0.75;
        std::vector<T> disp_j(ppc); disp_j[0] = 0.25; disp_j[1] = 0.75; disp_j[2] = 0.75; disp_j[3] = 0.25;
      
        int p = -1;
        for(int i = 0; i < NloopX; i++){
            for(int j = 0; j < NloopY; j++){
                for(int d = 0; d < ppc; d++){
                     p++;
                     T px = ((double)i+disp_i[d])*dx;
                     T py = ((double)j+disp_j[d])*dx;

                     particles.x[p](0) = px;
                     particles.x[p](1) = py;
		            particles.x0[p](0) = px;
		            particles.x0[p](1) = py;

                    } // end for d
                } // end for j
            } // end for i
        }

        for(int p = 0; p < Np; p++){
                  particles.v[p](0) = 0.;
                  particles.v[p](1) = 0.;

                  particles.theta[p] = slope_angle;
                  particles.rho[p] = material_density;
                  // particles.h[p] = 0.5;//sin(theta[p])*sqrt(py*py)+0.2;
                  particles.sigma_zz[p] = - 0.5*(particles.rho[p] /*- rho_water*/)*particles.h[p]*gravity*cos(particles.theta[p])*(1.-fluidization);

		  if (sim_type_int == 0) {
		     particles.sigma[p] = TM2::Zero();
		     particles.p0[p] = p0_Initial;
		     particles.epsilon_v_p_cumulative[p] = 1./xi * asinh(p0_Initial/K);
		  } else {
		    particles.sigma[p](0,0) = K0*particles.sigma_zz[p];
		    particles.sigma[p](1,1) = K0*particles.sigma_zz[p];
		    particles.sigma[p](0,1) = 0.5*(particles.sigma[p](0,0) - particles.sigma[p](1,1))*sin(phi);
		    particles.sigma[p](1,0) = 0.5*(particles.sigma[p](0,0) - particles.sigma[p](1,1))*sin(phi);
		  }

                  particles.mass[p] = particles.h[p] * particles.A[p]*particles.rho[p];
                  // particles.k[p] = xi*cohesion;
                  particles.surfacic_mass[p] = particles.h[p]*particles.rho[p];
                  particles.count[p] = 1.;


                  //For the release
		  particles.cohesion_WL[p] = Initial_cohesion_WL;
        //   if(particles.x[p](1) > 35.) {particles.cohesion_WL[p] = 1e6;}
        }



};

void Simulation::createDirectory(){

    if (mkdir(directory.c_str(), 0777) == -1)
        std::cerr << "Directory " << directory << " has already been created:  " << strerror(errno) << std::endl;
    else
        std::cout << "Directory " << directory << " was created now"<< std::endl;
    if (mkdir((directory + sim_name).c_str(), 0777) == -1)
        std::cerr << "Simulation " << sim_name << " has already been created: " << strerror(errno) << std::endl;
    else
        std::cout << "Simulation " << sim_name << " was created now" << std::endl;

}


void Simulation::simulate(){

    std::cout << "------------------------------------------------ " << std::endl;
    std::cout << "---------------- This is DAMPM ---------------- " << std::endl;
    std::cout << "------------------------------------------------ " << std::endl;

    createDirectory();

    // Write parameters to file for future reference
    // we need to adapt that in dampm with Lars' structure here
/*
    int emodel, pmodel;
    if (elastic_model == StvkWithHencky)
        emodel = 1;
    else if (elastic_model == NeoHookean)
        emodel = 2;
    else
        emodel = -1;
    if (plastic_model == VonMises)
        pmodel = 1;
    else if (plastic_model == NoPlasticity)
        pmodel = 2;
    else
        pmodel = -1;
*/
    if (sim_type_int == 0) {
        std::ofstream infoFile(directory + sim_name + "/info.txt");
        infoFile << end_frame               << "\n"  // 0
		 << frame_dt                << "\n"  // 1
		 << dx                      << "\n"  // 2
		 << rho                     << "\n"  // 3
		 << gravity                 << "\n"  // 4
		 << basal_frictional_angle  << "\n"  // 5
		 << Initial_cohesion_WL     << "\n"  // 6
		 << L1                      << "\n"  // 7
		 << L2                      << "\n"  // 8
		 << E                       << "\n"  // 9
		 << nu                      << "\n"  // 10
		 << Stiffness_WL            << "\n"  // 11
		 << velocity_PST            << "\n"  // 12
         << particles.theta[0]      << "\n"
         << p0_Initial              << "\n"
         << M                       << "\n"
         << beta                    << "\n"
         << flip_ratio              << "\n"
         << safety_factor           << "\n";
        infoFile.close();
    } else {
        std::ofstream infoFile(directory + sim_name + "/info.txt");
        infoFile << end_frame               << "\n"  // 0
                << frame_dt                << "\n"  // 1
                << dx                      << "\n"  // 2
                << rho                     << "\n"  // 3
                << gravity                 << "\n"  // 4
                << basal_frictional_angle  << "\n";  // 5
        infoFile.close();
    }


    // Total runtime of simulation
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Precomputations
    one_over_dx = 1.0 / dx;
    one_over_dx_square = one_over_dx * one_over_dx;

    // timing + saving
    time = 0;
    frame = 0;
    final_time = end_frame * frame_dt;
    saveParticle();
    // saveGridData();
    std::cout << "Particles saved" << std::endl;

    remesh();
    debug("Grid corner X", grid.xc);
    debug("Grid corner Y", grid.yc);

    // run
    while (frame < end_frame){
        std::cout << "Step: " << current_time_step << "    Time: " << time << std::endl;
        advanceStep();
        if (exit == 1)
            return;
        time += dt;
        current_time_step++;
        if( std::abs(time - frame_dt*(frame+1)) < 1e-15 ){
            frame++;
            std::cout << "Saving frame " << frame << std::endl;
            if(time > loading_time) {saveParticle();}
            // saveGridData();
        }
        if (std::abs(final_time-time) < 1e-15 || final_time < time){
            std::cout << "The simulation ended successfully at time = " << time << std::endl;
            break;
        }
    }
    saveParticle();
    // saveGridData();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Simulation took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " milliseconds" << std::endl;
    debug("Runtime P2G        = ", runtime_p2g     * 1000.0, " milliseconds");
    debug("Runtime G2P        = ", runtime_g2p     * 1000.0, " milliseconds");
    // debug("Runtime G2P1       = ", runtime_g2p_1     * 1000.0, " milliseconds");
    // debug("Runtime G2P2       = ", runtime_g2p_2     * 1000.0, " milliseconds");
    // debug("Runtime G2P3       = ", runtime_g2p_3     * 1000.0, " milliseconds");
    debug("Runtime comp for   = ", runtime_comp_for   * 1000.0, " milliseconds");
    debug("Runtime Elasticity = ", runtime_elasticity   * 1000.0, " milliseconds");
    debug("Runtime Plasticity = ", runtime_plasticity   * 1000.0, " milliseconds");
}

void Simulation::advanceStep(){

    updateDt();

    bool no_remesh = false;
    // ACTUAL (old) side lengths
    if (no_remesh == true) {
        grid.v.resize(Nx*Ny);    std::fill( grid.v.begin(),    grid.v.end(),    TV2::Zero() );
        grid.mass.resize(Nx*Ny); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );
	
    } else {
        remesh();
    }

    //

    P2G();
    compute_forces();
    compute_accel();
    G2P();
    if (parallelised)  {
      timer t_elasticity; t_elasticity.start();
      std::unique_ptr<std::vector<TM2>> stress_trial = elasticity_omp();
      t_elasticity.stop(); runtime_elasticity += t_elasticity.get_timing();
      // plasticity(stress_trial);
      //std::vector<TM2> stress_trial = elasticity_omp();
      timer t_plasticity; t_plasticity.start();
      plasticity_omp(stress_trial);
      t_plasticity.stop(); runtime_plasticity += t_plasticity.get_timing();
    } else {
      std::unique_ptr<std::vector<TM2>> stress_trial = elasticity();
      plasticity(stress_trial);
      //std::vector<TM2> stress_trial = elasticity();
      //plasticity(stress_trial);
    }

    // std::unique_ptr<std::vector<TM2>> stress_trial = elasticity();
    // plasticity(stress_trial);
    update_position();
    update_height();

    if (sim_type_int == 0) {
        // For the release (benchmark 1 & real field)
        // update_angle();

        // For benchmark 2 (PST)
        if (time > loading_time && time) {
            std::cout << "update_cohesion" << std::endl;
            update_cohesion();
        }
    }

}

void Simulation::compute_forces() {
  timer t_comp_for; t_comp_for.start();

    if (parallelised) {
        if (sim_type_int == 1) {
            compute_forces_flow_omp();
        } else {
            if (time <= loading_time) {
                compute_forces_release_loading_omp();
            } else {
                compute_forces_release_omp();
            }
        }
    } else {
        if (sim_type_int == 1) {
            compute_forces_flow();
        } else {
            if (time <= loading_time) {
                compute_forces_release_loading();
		//std::cout << "loading" << std::endl;
            } else {
                compute_forces_release();
		//std::cout << "loading finished" << std::endl;
            }
        }
    }
  t_comp_for.stop(); runtime_comp_for += t_comp_for.get_timing();
}

// void Simulation::update_cohesion(){
//   for (int p=0; p<Np; p++) {
//     TV2 xp;
//     xp(0) = particles.x0[p](0) - 250.;
//     xp(1) = particles.x0[p](1) - 250.;
//     if (xp.norm()  <= (time-loading_time) * velocity_PST){
//       particles.cohesion_WL[p] = 0.;
//     }
//   }
// }

void Simulation::update_cohesion(){
  for (int p=0; p<Np; p++) {
    TV2 xp;
    xp(0) = particles.x0[p](0);
    xp(1) = particles.x0[p](1);
    if (xp(1)  <= (time-loading_time) * velocity_PST){
      particles.cohesion_WL[p] = 0.;
    }
  }
}

void Simulation::update_angle(){
   for (int p=0; p<Np; p++) {
     particles.theta[p] = (M_PI/4.) + (time/final_time) * (M_PI/2.);
   }
}

void Simulation::update_height(){

  if (parallelised) {
    particles.h.resize(Np);    std::fill(particles.h.begin(),    particles.h.end(),    0.);
    particles.A.resize(Np);    std::fill(particles.A.begin(),    particles.A.end(),    0.);


#pragma omp parallel num_threads(n_threads)
    {
      std::vector<T> particles_h_local(Np, 0.);
      std::vector<T> particles_A_local(Np, 0.);
      //std::vector<T> particles_rho_local = particles.rho;

#pragma omp for
      for (int p = 0; p < Np; p++) {

            // rho update

        if (compressibility) {

	  T delta_epsilon_vol_plastic = (particles.I1_initial[p] - particles.I1[p])/K_3D;

	  particles.rho[p] /= (1. - delta_epsilon_vol_plastic);

        }

	// h update

	particles_h_local[p] = particles.surfacic_mass[p]/particles.rho[p];

	// A update

	particles_A_local[p] = particles.mass[p]/(particles_h_local[p]*particles.rho[p]);

      }
#pragma omp critical
      {
	for (int p = 0; p < Np; p++) {
	  particles.h[p]   += particles_h_local[p];
	  particles.A[p]   += particles_A_local[p];
	  //particles.rho[p] = particles_rho_local[p];
	}
      }
    }

  } else {
    for (int p = 0; p < Np; p++) {

      // rho update

      if (compressibility) {

	// T delta_epsilon_vol_plastic = (particles.I1_initial[p] - particles.I1[p])/K_3D;
	T delta_epsilon_vol_plastic = particles.d_epsilon_v_p[p];

	particles.rho[p] /= (1. - delta_epsilon_vol_plastic);

        }

      // h update

      particles.h[p] = particles.surfacic_mass[p]/particles.rho[p];

      // A update

      particles.A[p] = particles.mass[p]/(particles.h[p]*particles.rho[p]);

	 }
    }

}

void Simulation::compute_accel(){
    grid.a.resize(Nx*Ny); std::fill( grid.a.begin(), grid.a.end(), TV2::Zero() );
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            //NORMAL :
            if (grid.mass[ind(i,j,Ny)] > 1e-22) {
                grid.a[ind(i,j,Ny)] = grid.F[ind(i,j,Ny)]/grid.mass[ind(i,j,Ny)];
                grid.v[ind(i, j, Ny)] += dt*grid.a[ind(i, j, Ny)]; 
            }
        }
    }
}

void Simulation::P2G(){
    timer t_p2g; t_p2g.start();

    if (parallelised) {
        P2G_dampm_omp();
    } else {
        P2G_dampm();
    }
    // P2G_dampm_omp();
    t_p2g.stop(); runtime_p2g += t_p2g.get_timing();
} // end P2G



void Simulation::G2P(){
    timer t_g2p; t_g2p.start();

    if (parallelised)  {
        G2P_dampm_omp();
        // timer t_g2p_1; t_g2p_1.start();
        // G2P_dampm_omp_1();
        // t_g2p_1.stop();
        // timer t_g2p_2; t_g2p_2.start();
        // G2P_dampm_omp_2();
        // t_g2p_2.stop();
        // timer t_g2p_3; t_g2p_3.start();
        // G2P_dampm_omp_3();
        // t_g2p_3.stop();
        // runtime_g2p_1 += t_g2p_1.get_timing();
        // runtime_g2p_2 += t_g2p_2.get_timing();
        // runtime_g2p_3 += t_g2p_3.get_timing();
    } else {
        G2P_dampm();
    }
    // G2P_dampm_omp();

    t_g2p.stop(); runtime_g2p += t_g2p.get_timing();
}

void Simulation::update_position() {
    // std::fill(particles.D.begin(), particles.D.end(), TM2::Zero());

    for(int p=0; p<Np; p++) {
        particles.x[p] += dt*particles.pic[p];
        particles.v[p] = flip_ratio * ( particles.v[p] + particles.flip[p] ) + (1 - flip_ratio) * particles.pic[p];

        // // Update of the Inertia Matrix (APIC)
        // unsigned int i_base = std::floor((particles.x[p](0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        // unsigned int j_base = std::floor((particles.x[p](1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        // TV2 pos;
        // for(int i = i_base; i < i_base+4; i++){
        //     T xi = grid.x[i];

        //     for(int j = j_base; j < j_base+4; j++){
        //         T yi = grid.y[j];
        //         pos(0) = xi; pos(1) = yi;
        //         T weight = wip(particles.x[p](0), particles.x[p](1), xi, yi, one_over_dx);
        //         particles.D[p] += weight * (pos - particles.x[p]) * (pos - particles.x[p]).transpose();
        //     }
        // }
    }

}

void Simulation::updateDt(){
    // T dt_cfl = cfl * dx / max_speed;
    //  adaptive time step with max velocity
    // T max_speed = std::sqrt((particles.vx.array().square() + particles.vy.array().square()).maxCoeff());

    auto max_velocity_it = std::max_element( particles.v.begin(), particles.v.end(),
                                             []( const TV2 &v1, const TV2 &v2 )
                                             {
                                                 return v1.squaredNorm() < v2.squaredNorm();
                                             } );
    T max_speed = (*max_velocity_it).norm();
    // ///////////// DEBUG /////////////
    // debug("max_speed (iterator) = ", max_speed);
    // max_speed = -1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.v[p].norm() > max_speed){

    //         max_speed = particles.v[p].norm();
    //     }
    // }
    debug("max_speed      = ", max_speed);
    // /////////////////////////////////
    T dt_cfl = 0.;
    if (max_speed < 1e-10) {
        dt = dt_max;
        dt = std::min(dt, frame_dt*(frame+1) - time);
        dt = std::min(dt, final_time         - time);
      } else {
        dt_cfl = cfl * dx / max_speed;
        dt = std::min(dt_cfl, dt_max);
        dt = std::min(dt, frame_dt*(frame+1) - time);
        dt = std::min(dt, final_time         - time);
      }
    debug("               dt_cfl = ", dt_cfl);
    debug("               dt_max = ", dt_max);
    debug("               dt     = ", dt    );

    // if (dt > dt_cfl){
    //     debug("TIME STEP IS TOO BIG COMPARED TO CFL!!!");
    //     exit = 1;
    //     return;
    // }
    // if (dt > dt_max){
    //     debug("TIME STEP IS TOO BIG COMPARED TO ELASTIC WAVE SPEED!!!");
    //     exit = 1;
    //     return;
    // }
    // if (max_speed >= wave_speed){
    //     debug("DETECTED SPEED LARGER THAN ELASTIC WAVE SPEED!!!");
    //     exit = 1;
    //     return;
    // }
    // auto max_velocity_it = std::max_element( particles.v.begin(), particles.v.end(),
    //                                          []( const TV2 &v1, const TV2 &v2 )
    //                                          {
    //                                              return v1.squaredNorm() < v2.squaredNorm();
    //                                          } );
    // T max_speed = (*max_velocity_it).norm();
    // T dt = cfl * dx / max_speed;
    // // dt = 0.05;
    // dt = std::min(dt, dt_max);
    // dt = std::min(dt,frame_dt*(frame + 1)-time);
    // dt = std::min(dt,final_time - time);
} // end updateDt

void Simulation::remesh(){

    // ACTUAL min and max position of particles
    auto max_x_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV2 &x1, const TV2 &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T max_x = (*max_x_it)(0);
    auto max_y_it = std::max_element( particles.x.begin(), particles.x.end(),
                                             []( const TV2 &x1, const TV2 &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T max_y = (*max_y_it)(1);
    auto min_x_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV2 &x1, const TV2 &x2 )
                                             {
                                                 return x1(0) < x2(0);
                                             } );
    T min_x = (*min_x_it)(0);
    auto min_y_it = std::min_element( particles.x.begin(), particles.x.end(),
                                             []( const TV2 &x1, const TV2 &x2 )
                                             {
                                                 return x1(1) < x2(1);
                                             } );
    T min_y = (*min_y_it)(1);

    // ///////////// DEBUG /////////////
    // debug("max_x (iterator) = ", max_x);
    // debug("max_y (iterator) = ", max_y);
    // debug("min_x (iterator) = ", min_x);
    // debug("min_y (iterator) = ", min_y);
    //
    // max_x = -1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.x[p](0) > max_x){
    //         max_x = particles.x[p](0);
    //     }
    // }
    // debug("max_x (debug)    = ", max_x);
    //
    // max_y = -1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.x[p](1) > max_y){
    //         max_y = particles.x[p](1);
    //     }
    // }
    // debug("max_y (debug)    = ", max_y);
    //
    // min_x = 1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.x[p](0) < min_x){
    //         min_x = particles.x[p](0);
    //     }
    // }
    // debug("min_x (debug)    = ", min_x);
    //
    // min_y = 1e15;
    // for(int p=0; p<Np; p++){
    //     if (particles.x[p](1) < min_y){
    //         min_y = particles.x[p](1);
    //     }
    // }
    // debug("min_y (debug)    = ", min_y);
    // /////////////////////////////////

    // ACTUAL (old) side lengths
    T Lx = max_x - min_x;
    T Ly = max_y - min_y;

    // ACTUAL midpoints of particles
    T mid_x = min_x + Lx / 2.0;
    T mid_y = min_y + Ly / 2.0;

    // NEW number of grid-dx's per side length
    unsigned int safety_factor = 3;
    Nx = std::ceil(Lx * one_over_dx) + safety_factor;
    Ny = std::ceil(Ly * one_over_dx) + safety_factor;

    T low_x  = mid_x - Nx*dx/2.0;
    T low_y  = mid_y - Ny*dx/2.0;
    T high_x = mid_x + Nx*dx/2.0;
    T high_y = mid_y + Ny*dx/2.0;

    // NEW grid dimensions N = (L / dx + 1)
    Nx++;
    Ny++;

    debug("               grid   = (", Nx, ", ", Ny, ")");

    // Linspace x and y
    // Eigen:  LinSpaced(size, low, high) generates 'size' equally spaced values in the closed interval [low, high]
    grid.x = linspace(low_x, high_x, Nx);
    grid.y = linspace(low_y, high_y, Ny);

    grid.xc = grid.x[0]; // left bottom corner
    grid.yc = grid.y[0];

    grid.v.resize(Nx*Ny);    std::fill( grid.v.begin(),    grid.v.end(),    TV2::Zero() );
    // grid.flip.resize(Nx*Ny); std::fill( grid.flip.begin(), grid.flip.end(), TV2::Zero() );
    grid.mass.resize(Nx*Ny); std::fill( grid.mass.begin(), grid.mass.end(), 0.0 );
}

/*
void Simulation::moveObjects(T delta_t){
    for (InfinitePlate &obj : objects) {
        obj.move(delta_t);
    }
}
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// EXTRA FUNCTIONS //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// complex terrain tools:

TV2 Simulation::terrain_gradients(const int & p) { // computes 2D terrain gradients at particle p
    TV2 xp_old = particles.old_position[p];
    unsigned int i_base = std::floor((xp_old(0)-grid.xc)*one_over_dx); //we only want 2x2 square around particle
    unsigned int j_base = std::floor((xp_old(1)-grid.yc)*one_over_dx);

    // const unsigned int Nx = (grid.x).size();
    const unsigned int Ny = (grid.y).size();

    T dz_x = (grid.z[ind(i_base + 1,j_base,Ny)] + grid.z[ind(i_base + 1,j_base + 1,Ny)]
             - grid.z[ind(i_base,j_base,Ny)] - grid.z[ind(i_base,j_base +1,Ny)])/2./L1; // check L1 and L2 directions here <!>
    T dz_y = (grid.z[ind(i_base,j_base + 1,Ny)] + grid.z[ind(i_base + 1,j_base + 1,Ny)]
             - grid.z[ind(i_base + 1,j_base,Ny)] - grid.z[ind(i_base,j_base,Ny)])/2./L2;

    if (dz_x <= 1.e-6) { dz_x = 0.;}
    if (dz_y <= 1.e-6) { dz_y = 0.;}

    TV2 out;
    out << dz_x,
           dz_y;
    return out;

}
