#include "simulation.hpp"

void Simulation::G2P_dampm() {

    std::vector<TV2> v_L; v_L.resize(Nx*Ny);
    std::fill(v_L.begin(),    v_L.end(),    TV2::Zero());

    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            v_L[ind(i,j,Ny)] = grid.v[ind(i,j,Ny)] + dt*grid.a[ind(i,j,Ny)];
        }
    }

    grid.v.resize(Nx*Ny);    std::fill( grid.v.begin(),    grid.v.end(),    TV2::Zero() ); // velocity in 2D of nodes
    particles.pic.resize(Np); std::fill( particles.pic.begin(), particles.pic.end(), TV2::Zero() );
    particles.flip.resize(Np); std::fill( particles.flip.begin(), particles.flip.end(), TV2::Zero() );

    for (int p  = 0; p < Np; p++) {
        TV2 xp = particles.x[p];

        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        particles.a[p] = TV2::Zero();
        particles.old_position[p] = particles.x[p]; // copy old particles position

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x[i];

            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y[j];

        // compute vector wp (all interpolation weights for particule p)
            T weight = wip(xp(0), xp(1),xi,yi,one_over_dx);

            particles.pic[p]  += v_L[ind(i, j, Ny)]*weight;
            particles.flip[p] += dt*grid.a[ind(i, j, Ny)]*weight; 
            // particles.x[p] += dt*v_L[ind(i,j,Ny)]*weight;
            // particles.v[p] += dt*grid.a[ind(i,j,Ny)]*weight;
            // particles.a[p] += grid.a[ind(i,j,Ny)]*weight; // for info only

                }
            }

            for(int i = i_base; i < i_base+4; i++){
                T xi = grid.x[i];
                for(int j = j_base; j < j_base+4; j++){
                    T yi = grid.y[j];
                    T height_p = particles.h[p];
                    T m_p = particles.mass[p];
                    T weight = wip(xp(0), xp(1),xi,yi,one_over_dx); // use old positions
                    if (weight > 1e-22) {
                        grid.v[ind(i,j,Ny)] += m_p*(particles.v[p] + particles.flip[p])*weight; // updated particles.v at k+1
                    }
                }
            }
        // New particle velocity is a FLIP-PIC combination to be written

    } // end loop over particles


    for (int i = 0; i<Nx; i++) {
        for (int j = 0; j<Ny; j++) {
        ////// HERE MASS = 0
        if (grid.mass[ind(i,j,Ny)] > 0) {
            grid.v[ind(i,j,Ny)] /= grid.mass[ind(i,j,Ny)];
        } else {
            grid.v[ind(i,j,Ny)] = TV2::Zero();
        }
        if (grid.v[ind(i,j,Ny)] != grid.v[ind(i,j,Ny)]) {
            std::cout << "nan velocity g2p" << std::endl;
            std::cout <<  " mass : " << grid.mass[ind(i,j,Ny)] << std::endl;
        }
        }
    }
}
