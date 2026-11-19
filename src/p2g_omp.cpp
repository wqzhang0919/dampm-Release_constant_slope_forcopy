#include "simulation.hpp"
#include <omp.h>

void Simulation::P2G_dampm_omp(){

    // reset grid quantities
    grid.v.resize(Nx*Ny); std::fill( grid.v.begin(), grid.v.end(), TV2::Zero() );
    // grid.a.resize(Nx*Ny); std::fill( grid.a.begin(), grid.a.end(), TV2::Zero() );
    grid.mass.resize(Nx*Ny); std::fill( grid.mass.begin(), grid.mass.end(), 0.);


    #pragma omp parallel num_threads(n_threads)
    {
      std::vector<TV2> grid_v_local(Nx*Ny, TV2::Zero());
      //std::vector<TV2> grid_v_local(Nx*Ny);
      std::vector<T>   grid_mass_local(Nx*Ny, 0.);
      std::vector<TV2> grid_a_local(Nx*Ny, TV2::Zero());
	  
	  TM2 inv;

      #pragma omp for
      for(int p = 0; p < Np; p++){
		TV2 xp = particles.x[p];
		T m_p = particles.mass[p];
		// T height_p = particles.h[p];
        // T m_p = height_p*particles.rho[p]*particles.A[p];




	unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
	unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

	for(int i = i_base; i < i_base+4; i++){
	  T xi = grid.x[i];
	  for(int j = j_base; j < j_base+4; j++){
	    T yi = grid.y[j];
	    T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);
		//for APIC
		TV2 grad_weight = grad_wip(xp(0), xp(1), xi, yi, one_over_dx);

	    if (weight > 1e-25){

	      grid_mass_local[ind(i,j,Ny)]           += weight*m_p;  // w * h * A * rho
	    //   grid_mass_local[ind(i,j,Ny)]           += weight*m_p*particles.count[p];  // w * h * A * rho
	      grid_v_local[ind(i,j,Ny)]              += particles.v[p] * weight*m_p;//  + m_p*particles.B[p]*grad_weight; // OK
	    //   grid_a_local[ind(i,j,Ny)]              += particles.a[p] * weight; //*m_p;  // NOT mass weighted
	    }
	  } // end for j
	} // end for i
      } // end for p

      #pragma omp critical
      {
	for (int l = 0; l<Nx*Ny; l++){
	  grid.mass[l]      += grid_mass_local[l];
	  grid.v[l]         += grid_v_local[l];
	//   grid.a[l]         += grid_a_local[l];
	} // end for l
      } // end omp critical

    } // end omp parallel


    for (int i = 0; i<Nx; i++){
      for (int j = 0; j<Ny; j++) {
	T mi = grid.mass[ind(i,j,Ny)];
	if (mi > 1e-22) {
	  grid.v[ind(i,j,Ny)] /= mi;
	}
	else {
	  grid.v[ind(i,j,Ny)].setZero();
	}
      }// end for j
    } // end for i

} // end P2G_Optimized_Parallel
