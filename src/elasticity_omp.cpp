#include "simulation.hpp"
#include <omp.h>

std::unique_ptr< std::vector<TM2>> Simulation::elasticity_omp() {
    std::vector<TM2> strain_increment;
    strain_increment.resize(Np);    std::fill(strain_increment.begin(),    strain_increment.end(),    TM2::Zero());

    std::vector<TM2> stress_trial;
    stress_trial.resize(Np);    std::fill(stress_trial.begin(),    stress_trial.end(),    TM2::Zero());

    particles.surfacic_mass.resize(Np);    std::fill(particles.surfacic_mass.begin(),   particles.surfacic_mass.end(),    0.);

    #pragma omp parallel num_threads(n_threads)
    {

      //std::vector<TM2> strain_increment_local(Np, TM2::Zero());
      std::vector<T>   particles_surfacic_mass_local(Np, 0.);
      std::vector<TM2> stress_trial_local(Np, TM2::Zero());

      #pragma omp for
      for (int p = 0; p < Np; p++) {
        TV2 xp_old = particles.x[p];
        T height_p = particles.h[p];

        unsigned int i_base = std::floor((xp_old(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp_old(1)-grid.yc)*one_over_dx) - 1;

        for(int i = i_base; i < i_base + 4; i++) {
            T xi = grid.x[i];
            for(int j = j_base; j < j_base + 4; j++) {
	             T yi = grid.y[j];
	              if (grid.v[ind(i,j,Ny)] != grid.v[ind(i,j,Ny)]) {
		                std::cout << "nan velocity" << std::endl;
	              }
	              TM2 Gv_ip = grad_wip(xp_old(0),xp_old(1),xi,yi,one_over_dx)*grid.v[ind(i,j,Ny)].transpose(); // caution we need N{i,p} at time t^k and not t^{k+1} !!
	              strain_increment[p] += Gv_ip + Gv_ip.transpose();
                // strain_increment[p] *= particles.count[p];
            }
        }

        strain_increment[p] *= dt/2.0; // Delta epsilon

        T strain_trace = strain_increment[p].trace();

        if (strain_trace != strain_trace) {
	         strain_trace = 0.;
	          std::cout << " STRAIN TRACE NAN" << std::endl;
        }

	      // rho h update
	      particles_surfacic_mass_local[p] = particles.rho[p]*particles.h[p];

	      particles_surfacic_mass_local[p] /= (1. + strain_trace);


        // ELASTICITY

       stress_trial_local[p] = particles.sigma[p] + 2.*G*strain_increment[p] + (K-G)*strain_trace*TM2::Identity();
        // stress_trial_local[p] = particles.sigma[p] + E*strain_increment[p];
        //std::cout << "stress trial = " << stress_trial[p] << std::endl;
      }
    // return std::unique_ptr< std::vector<TM2> >(new stress_trial);

#pragma omp critical
      {
      	for (int p = 0; p < Np; p++) {
      	  particles.surfacic_mass[p] += particles_surfacic_mass_local[p];
      	  stress_trial[p]            += stress_trial_local[p];
      	}
      }
    } // end parallel

    return std::make_unique< std::vector<TM2> >(stress_trial);
    //return stress_trial;
}

// Remember:
// P = dPsidF              (first Piola-Kirchhoff stress tensor)
// tau = P * F.transpose() (Kirchhoff stress tensor)
