#include "simulation.hpp"

void Simulation::compute_forces_flow(){
    std::vector<TV2> b;
    b.resize(Np); std::fill(b.begin(),b.end(),TV2::Zero());
    
    // define that in Simulation?  

    // T dzxn; 
    // T dzyn;
    // T fluidization = 0;
    // T basal_frictional_angle = 25./180.*M_PI;

    grid.F.resize(Nx*Ny);    std::fill(grid.F.begin(),    grid.F.end(),    TV2::Zero());
    // NORMAL:
    for (int p=0; p<Np; p++) { // to be optimised here, store every weights and access via pointers + paralellize the loop on particles
        // first trial with pointers, needs to be reviewed
        // const TV2* xp = &particles.x[p];
        // const T* height_p = &particles.h[p];
        // const T m_p = *height_p*particles.rho[p]*particle_area;
        // std::cout << particles.rho[p] << std::endl;
        TV2 xp = particles.x[p];
        T height_p = particles.h[p];
        T A_p = particles.A[p];

        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        
        particles.sigma_xz[p] = 0.;//particles.rho[p]*(gravity)*cos(particles.theta[p])*height_p*(1.-fluidization)*tan(basal_frictional_angle);//0.;//particles.rho[p]*(-gravity)*cos(particles.theta[p])*height_p*(1.-fluidization)*tan(basal_frictional_angle);  // this has to be implemented when elevations.cpp is ready
        particles.sigma_yz[p] = particles.rho[p]*(-gravity)*cos(particles.theta[p])*height_p*(1.-fluidization)*tan(basal_frictional_angle);
        // particles.sigma_yz[p] = particles.rho[p]*( ((-gravity)*cos(particles.theta[p]) - abs(particles.v[p](1) * particles.v[p](1) * dthx[p] / yinc[p]) ) * (1.-fluidization)*tan(basal_frictional_angle) - abs(particles.v[p](1) * 3. * void_ratio * mu / height_p));
        // dthx(k) = datan(dzxn(njc(k))) -  datan(dzxn(jc(k))) // delta theta_x ? njc = new cell number where material point j is included
        //  acceleration centrifugal = vitesse angulaire  ^2 * r
        //      xinc(k) = dt*(vnx(n1)*s1(k) + vnx(n2)*s2(k) + vnx(n3)*s3(k) + vnx(n4)*s4(k)) // particle position increment = delta x p 

        particles.pressure[p] = 0.;//0.5*(particles.rho[p]- rho_water)*(-gravity*cos(particles.theta[p]))*height_p*fluidization + 0.5*rho_water*height_p*(-gravity*cos(particles.theta[p]));

        b[p](0) =  0.;//gravity*sin(particles.theta[p]);//0.;//-gravity*sin(particles.theta[p]); //sin(arctan(2D angle)) // this has to be implemented when elevations.cpp is ready
        b[p](1) = -gravity*sin(particles.theta[p]);
        
        for(int i = i_base; i < i_base + 4; i++){
            for(int j = j_base; j < j_base + 4; j++){
                    T xi = grid.x[i];
                    T yi = grid.y[j];
                    T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);
                    TV2 grad_weight = grad_wip(xp(0), xp(1), xi, yi, one_over_dx);

                grid.F[ind(i,j,Ny)](0) -= height_p*A_p*((particles.sigma[p](0,0) + particles.pressure[p])*grad_weight(0) + particles.sigma[p](0,1)*grad_weight(1) + particles.sigma_xz[p]*weight/height_p);
                grid.F[ind(i,j,Ny)](1) -= height_p*A_p*((particles.sigma[p](1,1) + particles.pressure[p])*grad_weight(1) + particles.sigma[p](0,1)*grad_weight(0) + particles.sigma_yz[p]*weight/height_p);
                
                grid.F[ind(i,j,Ny)](0) += height_p*A_p*particles.rho[p]*b[p](0)*weight;//+ 1./height_p)*weight;
                grid.F[ind(i,j,Ny)](1) += height_p*A_p*particles.rho[p]*b[p](1)*weight;

            }
        }
        
    }

}


void Simulation::compute_forces_release_loading(){
  std::vector<TV2> b;
  b.resize(Np); std::fill(b.begin(),b.end(),TV2::Zero());
    
    grid.F.resize(Nx*Ny);    std::fill(grid.F.begin(),    grid.F.end(),    TV2::Zero());
    T gravity_loading = gravity * (time / loading_time);


    // NORMAL:
    for (int p=0; p<Np; p++) { // to be optimised here, store every weights and access via pointers + paralellize the loop on particles
        // first trial with pointers, needs to be reviewed
        // const TV2* xp = &particles.x[p];
        // const T* height_p = &particles.h[p];
        // const T m_p = *height_p*particles.rho[p]*particle_area;
        // std::cout << particles.rho[p] << std::endl;
        TV2 xp = particles.x[p];
        T height_p = particles.h[p];
        T m_p = particles.mass[p];
	  
        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        TV2 displacement = particles.x[p] - particles.x0[p];
	T basal_stress_zz = particles.rho[p]*(gravity)*cos(particles.theta[p])*height_p;
	//T critical_displacement = (basal_stress_zz * tan(basal_frictional_angle) + Initial_cohesion_WL)/Stiffness_WL;
	T critical_displacement = Shear_strength/Stiffness_WL;
	
	TV2 direction_velocity = TV2::Zero() ;
	if (particles.v[p].norm() != 0) {
	  direction_velocity = particles.v[p]/(particles.v[p].norm());
	} else {
	  // std::cout << "No Velocity";
        }

	if ((abs(displacement(0)) <= critical_displacement) && (particles.cohesion_WL[p] > Initial_cohesion_WL/2.)) {
	  particles.sigma_xz[p] = Stiffness_WL * displacement(0);
	} else {
	  particles.cohesion_WL[p] = 0.;	  
	  particles.sigma_xz[p] = basal_stress_zz * tan(basal_frictional_angle) * (direction_velocity(0));
	  
	}

	if ((abs(displacement(1)) <= critical_displacement) && (particles.cohesion_WL[p] > Initial_cohesion_WL/2.)){
	  particles.sigma_yz[p] = Stiffness_WL * displacement(1);
	} else {
	  //std::cout << "NO COHESION";
	  particles.cohesion_WL[p] = 0.;
	  particles.sigma_yz[p] = basal_stress_zz * tan(basal_frictional_angle) * (direction_velocity(1));
	}
	
        b[p](0) =  0.;
        b[p](1) = -gravity_loading*sin(particles.theta[p]);

        for(int i = i_base; i < i_base+4; i++){
            T xi = grid.x[i];
            for(int j = j_base; j < j_base+4; j++){
                T yi = grid.y[j];
		T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);
                TV2 grad_weight = grad_wip(xp(0), xp(1), xi, yi, one_over_dx);
                if (grid.mass[ind(i,j,Ny)] > 0){
		  grid.F[ind(i,j,Ny)](0) -= height_p*particles.A[p]*((particles.sigma[p](0,0) + particles.pressure[p])*grad_weight(0) + particles.sigma[p](0,1)*grad_weight(1) + particles.sigma_xz[p]*weight/height_p);
		  grid.F[ind(i,j,Ny)](1) -= height_p*particles.A[p]*((particles.sigma[p](1,1) + particles.pressure[p])*grad_weight(1) + particles.sigma[p](0,1)*grad_weight(0) + particles.sigma_yz[p]*weight/height_p);
                
		  grid.F[ind(i,j,Ny)] += m_p*b[p]*weight;
		  // grid.F[ind(i,j,Ny)](1) += m_p*b[p](1)*weight;

		} // end if non-zero grid mass
             } // end for j
         } // end for i
    } // end for particles

}  

void Simulation::compute_forces_release(){
    std::vector<TV2> b;
    b.resize(Np); std::fill(b.begin(),b.end(),TV2::Zero());
    grid.F.resize(Nx*Ny);    std::fill(grid.F.begin(),    grid.F.end(),    TV2::Zero());


    // NORMAL:
    for (int p=0; p<Np; p++) {
        TV2 xp = particles.x[p];
        T height_p = particles.h[p];
        T m_p = particles.mass[p];

        unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
        unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines

        

	TV2 displacement = particles.x[p] - particles.x0[p];
	T basal_stress_zz = particles.rho[p]*(gravity)*cos(particles.theta[p])*height_p;
        //T critical_displacement = (basal_stress_zz * tan(basal_frictional_angle) + Initial_cohesion_WL)/Stiffness_WL;
	T critical_displacement = Shear_strength/Stiffness_WL;
	
	TV2 direction_velocity = TV2::Zero() ;
	if (particles.v[p].norm() != 0) {
	  direction_velocity = particles.v[p]/(particles.v[p].norm());
	} else {
	  std::cout << "No Velocity";
        }
	
	if ((abs(displacement(0)) <= critical_displacement) && (particles.cohesion_WL[p] > Initial_cohesion_WL/2.)) {
	  particles.sigma_xz[p] = Stiffness_WL * displacement(0);
    particles.sigma_xz[p] = 0.;
	} else {
	  particles.cohesion_WL[p] = 0.;	  
	  // particles.sigma_xz[p] = basal_stress_zz * tan(basal_frictional_angle) * (direction_velocity(0));
    particles.sigma_xz[p] = 0.;
	  
	}

	if ((abs(displacement(1)) <= critical_displacement) && (particles.cohesion_WL[p] > Initial_cohesion_WL/2.)){
	  particles.sigma_yz[p] = Stiffness_WL * displacement(1);
	} else {
	  particles.cohesion_WL[p] = 0.;
	  //particles.sigma_yz[p] = basal_stress_zz * tan(basal_frictional_angle) * (direction_velocity(1));
    particles.sigma_yz[p] = - basal_stress_zz * tan(basal_frictional_angle);
	}

	
        b[p](0) =  0.; //gravity*sin(particles.theta[p]);//0.;//-gravity*sin(particles.theta[p]); //sin(arctan(2D angle)) // this has to be implemented when elevations.cpp is ready
        b[p](1) = -gravity*sin(particles.theta[p]);

	                

         for(int i = i_base; i < i_base+4; i++) {
            for(int j = j_base; j < j_base+4; j++) {
                T xi = grid.x[i];
                    T yi = grid.y[j];
                    T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);
                    TV2 grad_weight = grad_wip(xp(0), xp(1), xi, yi, one_over_dx);
                if ((i >= i_base && i< i_base +4) && (j >= j_base && j < j_base + 4)){

                grid.F[ind(i,j,Ny)](0) -= height_p*particles.A[p]*((particles.sigma[p](0,0) + particles.pressure[p])*grad_weight(0) + particles.sigma[p](0,1)*grad_weight(1) + particles.sigma_xz[p]*weight/height_p);
                grid.F[ind(i,j,Ny)](1) -= height_p*particles.A[p]*((particles.sigma[p](1,1) + particles.pressure[p])*grad_weight(1) + particles.sigma[p](0,1)*grad_weight(0) + particles.sigma_yz[p]*weight/height_p);
                
                grid.F[ind(i,j,Ny)] += m_p*b[p]*weight;//+ 1./height_p)*weight;
		// grid.F[ind(i,j,Ny)](1) += m_p*b[p](1)*weight;

                }
            }
        }
        
    }


}
