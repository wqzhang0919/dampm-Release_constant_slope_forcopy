#include "simulation.hpp"
#include <omp.h>

void Simulation::G2P_dampm_omp() {
  particles.old_position = particles.x;
  particles.pic.resize(Np); std::fill( particles.pic.begin(), particles.pic.end(), TV2::Zero() );
  particles.flip.resize(Np); std::fill( particles.flip.begin(), particles.flip.end(), TV2::Zero() );

  #pragma omp parallel num_threads(n_threads)
  {
    std::vector<TV2> particles_x_local(Np, TV2::Zero() );
    std::vector<TV2> particles_pic_local(Np, TV2::Zero() );
    std::vector<TV2> particles_flip_local(Np, TV2::Zero() );
    // std::vector<TV2> particles_a_local(Np, TV2::Zero() );

    // std::vector<TV2> particles_v_local(Np);
    // std::vector<TV2> particles_x_local(Np);
    // std::vector<TV2> particles_a_local(Np);

      #pragma omp for
      for(int p = 0; p < Np; p++){
          TV2 xp = particles.x[p];

    unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
    unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both 	quadratic and cubic splines

    for(int i = i_base; i < i_base+4; i++){
      T xi = grid.x[i];
      for(int j = j_base; j < j_base+4; j++){
        T yi = grid.y[j];
        T weight = wip(xp(0), xp(1),xi,yi,one_over_dx);

        // particles_x_local[p] += dt*v_L[ind(i,j,Ny)]*weight;
        // particles_x_local[p] += dt*(grid.v[ind(i,j,Ny)] + dt*grid.a[ind(i,j,Ny)])*weight;
        // particles_pic_local[p] += (grid.v[ind(i,j,Ny)] + dt*grid.a[ind(i,j,Ny)])*weight;
        particles_pic_local[p] += grid.v[ind(i,j,Ny)]*weight;
        particles_flip_local[p] += dt*grid.a[ind(i, j, Ny)]*weight;
        // particles_a_local[p] += grid.a[ind(i,j,Ny)]*weight; // for info only

      } // end loop j
    } // end loop i

  } // end loop p

      #pragma omp critical
      {
        for (int p = 0; p<Np; p++) {
          // particles.x[p]            += particles_x_local[p];
          particles.pic[p] += particles_pic_local[p];
          particles.flip[p] += particles_flip_local[p];
        }
      } // end omp critical

  } //end opm parallel

    //   //Interpolation back of B matrix (APIC)
    // particles.B.resize(Np); std::fill( particles.B.begin(), particles.B.end(), TM2::Zero() );
    //  #pragma omp parallel num_threads(n_threads)
    // {
    //  std::vector<TM2> particles_B_loc(Np, TM2::Zero() );
    //  #pragma omp for
    //  for (int p = 0; p<Np; p++) {

    //    TV2 xp = particles.x[p];

    //    unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
    //    unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
    //     TV2 pos;
    //    for(int i = i_base; i < i_base+4; i++){
	  //       T xi = grid.x[i];

	  //       for(int j = j_base; j < j_base+4; j++){
	  //         T yi = grid.y[j];
    //         pos(0) = xi; pos(1) = yi;
	  //         T weight = wip(xp(0), xp(1),xi,yi,one_over_dx); // use old positions
	  //         if (weight > 1e-25) {
    //           particles_B_loc[p] += weight * grid.v[ind(i, j, Ny)] * (pos - particles.x[p]).transpose();
	  //         }
	  //       }
    //    }
    //  }

    //  #pragma omp critical
    //  {

    //    for (int p = 0; p<Np; p++) {
    //      particles.B[p] += particles_B_loc[p];
    //    }
    //  }
    // }


  // grid.v.resize(Np); std::fill( grid.v.begin(), grid.v.end(), TV2::Zero() );
  //  #pragma omp parallel num_threads(n_threads)
  //   {
  //    std::vector<TV2> grid_v_loc(Nx*Ny, TV2::Zero() );
  //    #pragma omp for
  //    for (int p = 0; p<Np; p++) {

  //      TV2 xp = particles.x[p];

  //      unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
  //      unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

  //      for(int i = i_base; i < i_base+4; i++){
	//         T xi = grid.x[i];

	//         for(int j = j_base; j < j_base+4; j++){
	//           T yi = grid.y[j];
	//           T m_p = particles.mass[p];
	//           T weight = wip(xp(0), xp(1),xi,yi,one_over_dx); // use old positions
  //           TV2 grad_weight = grad_wip(xp(0), xp(1),xi,yi,one_over_dx);
	//           if (weight > 1e-25) {
	//             // grid_v_loc[ind(i,j,Ny)] += m_p * (particles.v[p] + particles.flip[p]) * weight; // updated particles.v at k+1
	//             grid_v_loc[ind(i,j,Ny)] += m_p * (particles.v[p] + particles.flip[p])* weight;// + m_p*particles.B[p]*grad_weight;
  //             // grid_v_loc[ind(i,j,Ny)] += m_p * (particles.v[p] + particles.flip[p]) * weight*particles.count[p]; // updated particles.v at k+1
	//           }
	//         }
  //      }
  //    }

  //    #pragma omp critical
  //    {

  //      for (int l = 0; l<Nx*Ny; l ++) {
	//         grid.v[l] += grid_v_loc[l];
  //      }

  //    }
  //   }
  //   // end omp parallel


  //   for (int i = 0; i<Nx; i++) {
  //     for (int j = 0; j<Ny; j++) {
	//       ////// HERE MASS = 0
	//       if (grid.mass[ind(i,j,Ny)] > 0) {
	//         grid.v[ind(i,j,Ny)] /= grid.mass[ind(i,j,Ny)];
	//       } else {
	//         grid.v[ind(i,j,Ny)] = TV2::Zero();
	//       }
	//       if (grid.v[ind(i,j,Ny)] != grid.v[ind(i,j,Ny)]) {
	//         std::cout << "nan velocity g2p" << std::endl;
	//         std::cout <<  " mass : " << grid.mass[ind(i,j,Ny)] << std::endl;
	//       }
  //     }
  //   }


} //END FUNCTION

//
// void Simulation::G2P_dampm_omp(){
//
//     std::vector<TV2> v_L; v_L.resize(Nx*Ny);
//     std::fill(v_L.begin(),    v_L.end(),    TV2::Zero());
//
//     #pragma omp parallel num_threads(n_threads)
//     {
//       std::vector<TV2> v_L_local(Np, TV2::Zero() );
//
//       #pragma omp for
//       for (int i = 0; i < Nx; i++){
//         for (int j = 0; j < Ny; j++){
//             v_L_local[ind(i,j,Ny)] = grid.v[ind(i,j,Ny)] + dt*grid.a[ind(i,j,Ny)];
//         }
//       }
//       #pragma omp critical
//       {
//         for (int l = 0; l<Nx*Ny; l ++) {
//            v_L[l] += v_L_local[l];
//         }
//       }
//     } // end omp parallel
//
//     particles.old_position = particles.x;
//     // for(int p = 0; p < Np; p++){
//     //   particles.old_position[p] = particles.x[p];
//     // }
//
//     grid.v.resize(Nx*Ny);
//     std::fill( grid.v.begin(),    grid.v.end(),    TV2::Zero() );
//     std::fill( particles.a.begin(), particles.a.end(), TV2::Zero() );
//
//
//     #pragma omp parallel num_threads(n_threads)
//     {
//       std::vector<TV2> particles_v_local(Np, TV2::Zero() );
//       std::vector<TV2> particles_x_local(Np, TV2::Zero() );
//       std::vector<TV2> particles_a_local(Np, TV2::Zero() );
//
//       // std::vector<TV2> particles_v_local(Np);
//       // std::vector<TV2> particles_x_local(Np);
//       // std::vector<TV2> particles_a_local(Np);
//
//         #pragma omp for
//         for(int p = 0; p < Np; p++){
//             TV2 xp = particles.x[p];
//
// 	    unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
// 	    unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both 	quadratic and cubic splines
//
// 	    for(int i = i_base; i < i_base+4; i++){
// 	      T xi = grid.x[i];
// 	      for(int j = j_base; j < j_base+4; j++){
//       		T yi = grid.y[j];
// 		T weight = wip(xp(0), xp(1),xi,yi,one_over_dx);
//
// 		particles_x_local[p] += dt*v_L[ind(i,j,Ny)]*weight;
// 		particles_v_local[p] += dt*grid.a[ind(i,j,Ny)]*weight;
// 		particles_a_local[p] += grid.a[ind(i,j,Ny)]*weight; // for info only
//
// 	      } // end loop j
//             } // end loop i
//
//         } // end loop p
//
//         #pragma omp critical
//         {
// 	  for (int p = 0; p<Np; p++) {
// 	    particles.x[p]            += particles_x_local[p];
// 	    particles.v[p]            += particles_v_local[p];
// 	    particles.a[p]            += particles_a_local[p];
// 	  }
//         } // end omp critical
//
//     } //end opm parallel
//
//
//
//     #pragma omp parallel num_threads(n_threads)
//     {
//      std::vector<TV2> grid_v_loc(Nx*Ny, TV2::Zero() );
//      #pragma omp for
//      for (int p = 0; p<Np; p++) {
//
//        TV2 xp = particles.old_position[p];
//
//        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
//        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
//
//        for(int i = i_base; i < i_base+4; i++){
// 	 T xi = grid.x[i];
//
// 	 for(int j = j_base; j < j_base+4; j++){
// 	   T yi = grid.y[j];
// 	   T m_p = particles.mass[p];
// 	   T weight = wip(xp(0), xp(1),xi,yi,one_over_dx); // use old positions
// 	   if (weight > 1e-25) {
// 	     grid_v_loc[ind(i,j,Ny)] += m_p * particles.v[p] * weight; // updated particles.v at k+1
// 	   }
// 	 }
//        }
//      }
//
//      #pragma omp critical
//      {
//
//        for (int l = 0; l<Nx*Ny; l ++) {
// 	        grid.v[l] += grid_v_loc[l];
//        }
//
//      }
//     }
//     // end omp parallel
//
//
//     for (int i = 0; i<Nx; i++) {
//       for (int j = 0; j<Ny; j++) {
// 	////// HERE MASS = 0
// 	if (grid.mass[ind(i,j,Ny)] > 0) {
// 	  grid.v[ind(i,j,Ny)] /= grid.mass[ind(i,j,Ny)];
// 	} else {
// 	  grid.v[ind(i,j,Ny)] = TV2::Zero();
// 	}
// 	if (grid.v[ind(i,j,Ny)] != grid.v[ind(i,j,Ny)]) {
// 	  std::cout << "nan velocity g2p" << std::endl;
// 	  std::cout <<  " mass : " << grid.mass[ind(i,j,Ny)] << std::endl;
// 	}
//       }
//     }
//
// } // end G2P_Optimized_Parallel
