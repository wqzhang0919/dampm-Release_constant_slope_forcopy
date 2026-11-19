#include "simulation.hpp"

void Simulation::P2G_dampm() {

    // reset grid quantities
        grid.v.resize(Nx*Ny); std::fill( grid.v.begin(), grid.v.end(), TV2::Zero() );
        grid.a.resize(Nx*Ny); std::fill( grid.a.begin(), grid.a.end(), TV2::Zero() );
        grid.mass.resize(Nx*Ny); std::fill( grid.mass.begin(), grid.mass.end(), 0.);

    // interpolation
    for(int p = 0; p < Np; p++){
        TV2 xp = particles.x[p];
        T height_p = particles.h[p];
        T m_p = particles.mass[p]; 
        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x[i];

            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y[j];

                T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);
                if (weight > 1e-22){
                
                    grid.mass[ind(i,j,Ny)]           += weight*m_p;  // w * h * A * rho
                    grid.v[ind(i,j,Ny)]              += particles.v[p] * weight*m_p; // OK 
                    grid.a[ind(i,j,Ny)]              += particles.a[p] * weight; //*m_p;  // NOT mass weighted
                }

            } // end for j
        } // end for i
    } // end for p

    for (int i = 0; i<Nx; i++){
        for (int j = 0; j<Ny; j++) {
             T mi = grid.mass[ind(i,j,Ny)];
            if (mi > 1e-22) {
                grid.v[ind(i,j,Ny)] /= mi;
                }
            else {
                grid.v[ind(i,j,Ny)].setZero(); 
            }

        }
       
    }

}
