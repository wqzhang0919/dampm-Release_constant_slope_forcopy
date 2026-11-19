#include "simulation.hpp"
#include "tools.hpp"
using namespace std;


void Simulation::saveParticle(string extra) {
  if (saving_type == "ALL") {
    saveParticleData();
  } else if (saving_type == "Light") {
    saveParticleData_light();
  } else if (saving_type == "PST") {
    saveParticleData_PST();
  } else if (saving_type == "CrossSlope") {
    saveParticleData_2D();
  }

}

void Simulation::saveParticleData_PST(string extra){

  ofstream outFile(directory + sim_name + "/out_part_frame_" + extra + to_string(frame) + ".csv");

  outFile    << "y"           << ","   // 0
	     << "vx"          << ","   // 1
	     << "vy"          << ","   // 2
	     << "ux"          << ","   // 3
  	  << "uy"          << ","   // 4
	     << "cohesion_WL" << ","   // 5
	     << "sigma_yz"    << ","   // 6
	     << "sigma_xz"    << ","   // 7
	     << "sigma_yy"    << ","   // 8
	     << "sigma_xx"    << ","   // 9
	     << "sigma_xy"    << ","   // 10
	     << "h"           << "\n"; // 11

  unsigned int NloopX = round(L1/dx);
  unsigned int NloopY = round(L2/dx);
  unsigned int ppc = 4.;

  cout << NloopY << endl;

  for (int j=0; j<NloopY; j++) {

    int p = ppc*j;
    TV2 displacement = particles.x[p] - particles.x0[p];

    T mean_y = 0.;
    T mean_vx = 0.;
    T mean_vy = 0.;
    T mean_ux = 0.;
    T mean_uy = 0.;
    T mean_Cohesion = 0.;
    T mean_sigma_yz = 0.;
    T mean_sigma_xz = 0.;
    T mean_sigma_yy = 0.;
    T mean_sigma_xx = 0.;
    T mean_sigma_xy = 0.;
    T mean_h = 0.;

    // T mean_y = particles.x[p](1);
    // T mean_vx = particles.v[p](0);
    // T mean_vy = particles.v[p](1);
    // T mean_ux = displacement(0);
    // T mean_uy = displacement(1);
    // T mean_Cohesion = particles.cohesion_WL[p];
    // T mean_sigma_yz = particles.sigma_yz[p];
    // T mean_sigma_xz = particles.sigma_xz[p];
    // T mean_sigma_yy = particles.sigma[p](1, 1);
    // T mean_sigma_xx = particles.sigma[p](0, 0);
    // T mean_sigma_xy = particles.sigma[p](0, 1);

    for (int i=0; i<NloopX; i++) {
      for (int d=0; d<ppc; d++) {
	p = ppc*j + i*ppc*NloopY + d;

	TV2 displacement = particles.x[p] - particles.x0[p];

        mean_y        += particles.x[p](1);
	mean_vx       += particles.v[p](0);
	mean_vy       += particles.v[p](1);
	mean_ux       += displacement(0);
	mean_uy       += displacement(1);
	mean_Cohesion += particles.cohesion_WL[p];
	mean_sigma_yz += particles.sigma_yz[p];
	mean_sigma_xz += particles.sigma_xz[p];
	mean_sigma_yy += particles.sigma[p](1, 1);
	mean_sigma_xx += particles.sigma[p](0, 0);
	mean_sigma_xy += particles.sigma[p](0, 1);
	mean_h        += particles.h[p];

      }
    }
    mean_y        /= ppc*NloopX;
    mean_vx       /= ppc*NloopX;
    mean_vy       /= ppc*NloopX;
    mean_ux       /= ppc*NloopX;
    mean_uy       /= ppc*NloopX;
    mean_Cohesion /= ppc*NloopX;
    mean_sigma_yz /= ppc*NloopX;
    mean_sigma_xz /= ppc*NloopX;
    mean_sigma_yy /= ppc*NloopX;
    mean_sigma_xx /= ppc*NloopX;
    mean_sigma_xy /= ppc*NloopX;
    mean_h        /= ppc*NloopX;

    outFile << mean_y          << ","   // 0
	    << mean_vx         << ","   // 1
	    << mean_vy         << ","   // 2
	    << mean_ux         << ","   // 3
	    << mean_uy         << ","   // 4
	    << mean_Cohesion   << ","   // 5
	    << mean_sigma_yz   << ","   // 6
	    << mean_sigma_xz   << ","   // 7
	    << mean_sigma_yy   << ","   // 8
	    << mean_sigma_xx   << ","   // 9
	    << mean_sigma_xy   << ","   // 10
	    << mean_h          << "\n"; // 11

  }

}

void Simulation::saveParticleData(string extra){
    ofstream outFile(directory + sim_name + "/out_part_frame_" + extra + to_string(frame) + ".csv");

    outFile     << "x"           << ","   // 0
                << "y"           << ","   // 1
		            << "vx"          << ","   // 2
		            << "vy"          << ","   // 3
                // << "ux"          << ","
                // << "uy"          << ","
                << "h"           << ","
      //        << "dt"          << ","
                << "theta"       << ","
                << "cohesion_WL" << ","
                << "sigma_xz"    << ","
                << "sigma_yz"    << ","
                << "sigma_yy"    << ","
                << "sigma_xx"    << ","
                << "sigma_xy"    << ","
                // << "x_proj" << ","
                // << "y_proj" << ","
                << "p0"     << "\n";
                // << "area"        << ","
                // << "ax"          << ","
                // << "ay"          << ","
                // // << "I1_initial"  << ","
                // // << "J2_initial"  << ","
                // // << "k_initial"   << ","
                // // << "k_final"     << ","
                // // << "I1"          << ","
                // // << "J2"          << ","
                // << "f"           << ","
                // << "RMA_iter"    << ',' //<< "\n";   // 4
                // << "gradx" << ","
                // << "grady" << "\n";

    for(int p = 0; p < Np; p++){
        TV2 displacement = particles.x[p] - particles.x0[p];
        // outFile << setprecision(numeric_limits<long double>::digits10 + 1)
	outFile       << particles.x[p](0)          << ","   // 0
                << particles.x[p](1)          << ","   // 1
	              << particles.v[p](0)          << ","   // 3
	              << particles.v[p](1)          << ","
	              // << displacement(0)            << ","   // 4
	  	          // << displacement(1)            << ","   // 5
                << particles.h[p]             << ","
	  //    << dt                         << ","
	              << particles.theta[p]         << ","   // 8
                << particles.cohesion_WL[p]   << ","   // 9
	              << particles.sigma_xz[p]      << ","   // 10
                << particles.sigma_yz[p]      << ","   // 11
                << particles.sigma[p](1,1)    << ","   // 12
                << particles.sigma[p](0,0)    << ","   // 13
                << particles.sigma[p](0,1)    << ","   // 14
                // << particles.x_proj[p](0) << ","
                // << particles.x_proj[p](1) << ","
                <<particles.p0[p] << "\n";
                // << particles.A[p]             << ","   // 15
                // << particles.a[p](0)          << ","   // 3
                // << particles.a[p](1)          << ","
	        // // << particles.I1_initial[p]    << ","
                // // << particles.J2_initial[p]    << ","
                // // << particles.k_initial[p]     << ","
                // // << particles.k_final[p]       << ","
                // // << particles.I1[p]            << ","
                // // << particles.J2[p]            << ","
                // << particles.f[p]             << ","
                // << particles.count[p]         << ',' // << "\n";  // 4
                // << particle_gradients(p, false)(0) << ","
                // << particle_gradients(p, false)(1) <<  "\n";
    }
}

void Simulation::saveParticleData_light(string extra){
    ofstream outFile(directory + sim_name + "/out_part_frame_" + extra + to_string(frame) + ".csv");

    outFile     << "x"           << ","   // 0
                << "y"           << ","   // 1
		            << "vx"          << ","   // 2
		            << "vy"          << ","   // 3
      //        << "ux"          << ","
      //        << "uy"          << ","
                << "h"           << ","
      //        << "dt"          << ","
      //        << "theta"       << ","
                << "cohesion_WL" << ","
                // << "sigma_xz"    << ","
                // << "sigma_yz"    << ","
                << "tau_0" << ","
                << "sigma_yy"    << ","
                << "sigma_xx"    << ","
                << "sigma_xy"    << ","
                << "p0"     << "\n";
                // << "area"        << ","
                // << "ax"          << ","
                // << "ay"          << ","
                // // << "I1_initial"  << ","
                // // << "J2_initial"  << ","
                // // << "k_initial"   << ","
                // // << "k_final"     << ","
                // // << "I1"          << ","
                // // << "J2"          << ","
                // << "f"           << ","
                // << "RMA_iter"    << ',' //<< "\n";   // 4
                // << "gradx" << ","
                // << "grady" << "\n";

    for(int p = 0; p < Np; p++){
        // TV2 displacement = particles.x[p] - particles.x0[p];
        T tau_0 = sqrt(particles.sigma_xz[p]*particles.sigma_xz[p] + particles.sigma_yz[p]*particles.sigma_yz[p]);
        // outFile << setprecision(numeric_limits<long double>::digits10 + 1)
	outFile << particles.x[p](0)          << ","   // 0
                << particles.x[p](1)          << ","   // 1
	        << particles.v[p](0)          << ","   // 3
	        << particles.v[p](1)          << ","
	  //    << displacement(0)            << ","   // 4
	  // 	<< displacement(1)            << ","   // 5
                << particles.h[p]             << ","
	  //    << dt                         << ","
	  //    << particles.theta[p]         << ","   // 8
                << particles.cohesion_WL[p]   << ","   // 9
	  //    << particles.sigma_xz[p]      << ","   // 10
                // << particles.sigma_yz[p]      << ","   // 11
                << tau_0 << ","
                << particles.sigma[p](1,1)    << ","   // 12
                << particles.sigma[p](0,0)    << ","   // 13
                << particles.sigma[p](0,1)    << ","   // 14
                <<particles.p0[p] << "\n";
                // << particles.A[p]             << ","   // 15
                // << particles.a[p](0)          << ","   // 3
                // << particles.a[p](1)          << ","
	        // // << particles.I1_initial[p]    << ","
                // // << particles.J2_initial[p]    << ","
                // // << particles.k_initial[p]     << ","
                // // << particles.k_final[p]       << ","
                // // << particles.I1[p]            << ","
                // // << particles.J2[p]            << ","
                // << particles.f[p]             << ","
                // << particles.count[p]         << ',' // << "\n";  // 4
                // << particle_gradients(p, false)(0) << ","
                // << particle_gradients(p, false)(1) <<  "\n";
    }
}

void Simulation::saveGridData(std::string extra){
    std::ofstream outFile(directory + sim_name + "/out_grid_frame_" + extra + std::to_string(frame) + ".csv");
    outFile         << "x"       << "," << "y"       << "," << "vx"         << "," << "vy"         << ","  << "mass"         << "," << "Fx" << "," << "Fy"  << "," << "ax" << "," << "ay" << "\n";
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            outFile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << grid.x[i] << "," << grid.y[j] << "," << grid.v[ind(i,j,Ny)](0) << "," << grid.v[ind(i,j,Ny)](1) << "," << grid.mass[ind(i,j,Ny)] << ","  <<  grid.F[ind(i,j,Ny)](0) << "," << grid.F[ind(i,j,Ny)](1) << "," << grid.a[ind(i,j,Ny)](0) << "," << grid.a[ind(i,j,Ny)](1) <<  "\n";
        }
    }
}




void Simulation::saveParticleData_2D(std::string extra){
  std::ofstream outFile(directory + sim_name + "/out_part_vertical_frame_" + extra + std::to_string(frame) + ".csv");
  
  outFile     << "x"           << ","   // 0
	      << "y"           << ","   // 1
	      << "vx"          << ","   // 2
	      << "vy"          << ","   // 3
	      << "h"           << ","
	      << "dt"          << ","
	      << "cohesion_WL" << ","
	      << "sigma_yz"    << ","
	      << "sigma_xz"    << ","
	      << "sigma_yy"    << ","
	      << "sigma_xx"    << ","
	      << "sigma_xy"    << ","
	      << "area"        << "\n";

  unsigned int NloopX = std::round(L1/dx);
  unsigned int NloopY = std::round(L2/dx);
  unsigned int ppc = 4.;
  int p;
  
  for(int j = 0; j < 2*NloopY; j++){
    p = 2*j + std::round(NloopX/2.)*ppc*NloopY;
    
    //outFile << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
    outFile << particles.x[p](0)          << ","   // 0
	    << particles.x[p](1)          << ","   // 1
	    << particles.v[p](0)          << ","   // 3
	    << particles.v[p](1)          << ","
	    << particles.h[p]             << ","
	    << dt                         << ","
	    << particles.cohesion_WL[p]   << ","   // 9
	    << particles.sigma_yz[p]      << ","   // 10
	    << particles.sigma_xz[p]      << ","   // 11
	    << particles.sigma[p](1,1)    << ","   // 12
	    << particles.sigma[p](0,0)    << ","   // 13
	    << particles.sigma[p](0,1)    << ","   // 14
	    << particles.A[p]             << "\n";   // 15
  }


    std::ofstream outFile1(directory + sim_name + "/out_part_horizontal_frame_" + extra + std::to_string(frame) + ".csv");
  
  outFile1     << "x"           << ","   // 0
	      << "y"           << ","   // 1
	      << "vx"          << ","   // 2
	      << "vy"          << ","   // 3
	      << "h"           << ","
	      << "dt"          << ","
	      << "cohesion_WL" << ","
	      << "sigma_yz"    << ","
	      << "sigma_xz"    << ","
	      << "sigma_yy"    << ","
	      << "sigma_xx"    << ","
	      << "sigma_xy"    << ","
	      << "area"        << "\n";
  
  for(int i = 0; i < NloopX; i++){
    p = i*ppc*NloopY + std::round(NloopY/2.)*ppc;
    
    //outFile << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
    outFile1 << particles.x[p](0)          << ","   // 0
	    << particles.x[p](1)          << ","   // 1
	    << particles.v[p](0)          << ","   // 3
	    << particles.v[p](1)          << ","
	    << particles.h[p]             << ","
	    << dt                         << ","
	    << particles.cohesion_WL[p]   << ","   // 9
	    << particles.sigma_yz[p]      << ","   // 10
	    << particles.sigma_xz[p]      << ","   // 11
	    << particles.sigma[p](1,1)    << ","   // 12
	    << particles.sigma[p](0,0)    << ","   // 13
	    << particles.sigma[p](0,1)    << ","   // 14
	    << particles.A[p]             << "\n";   // 15

    p = i*ppc*NloopY + std::round(NloopY/2.)*ppc+ 3;
    outFile1 << particles.x[p](0)          << ","   // 0
	    << particles.x[p](1)          << ","   // 1
	    << particles.v[p](0)          << ","   // 3
	    << particles.v[p](1)          << ","
	    << particles.h[p]             << ","
	    << dt                         << ","
	    << particles.cohesion_WL[p]   << ","   // 9
	    << particles.sigma_yz[p]      << ","   // 10
	    << particles.sigma_xz[p]      << ","   // 11
	    << particles.sigma[p](1,1)    << ","   // 12
	    << particles.sigma[p](0,0)    << ","   // 13
	    << particles.sigma[p](0,1)    << ","   // 14
	    << particles.A[p]             << "\n";   // 15
  }
}
