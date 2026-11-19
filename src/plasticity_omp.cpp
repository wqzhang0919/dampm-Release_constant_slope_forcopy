#include "simulation.hpp"
#include <omp.h>

void Simulation::plasticity_omp(std::unique_ptr< std::vector<TM2>> & stress_trial) {


  if (plasticity_model_int == 1) { // Drucker_prager
    auto f_dp  = [](T &J2, T &p_I1, T &k, T eta) { return sqrt(J2) + eta*p_I1 - k; }; // see De Souza

#pragma omp parallel num_threads(n_threads)
    {
#pragma omp for
      for (int p = 0; p < Np; p++) {
        // PLASTICITY : DRUCKER-PRAGER

        // zz component of sigma

        particles.sigma_zz[p] = - 0.5*particles.rho[p]*particles.h[p]*gravity*cos(particles.theta[p]);

        // initial hardening constant
        particles.k_initial[p] = particles.k[p];

        // first invariant of Cauchy stress tensor
        T I1 = (*stress_trial)[p].trace() + particles.sigma_zz[p];

        T p_I1 = I1/2.;
	
        std::cout << "I1  = " << I1 << std::endl;

        // deviatoric Cauchy stress tensor
        TM2 s = (*stress_trial)[p] - p_I1*TM2::Identity();
        std::cout << "s = " << s << std::endl;

        T s_zz = particles.sigma_zz[p] - p_I1;
        // TA2 s_array = s;

        // second invariant of deviatoric stress tensor
        T J2 = s(0,1)*s(0,1) - s(0,0)*s(1,1) - s(0,0)*s_zz - s(1,1)*s_zz + 1.e-15;

        particles.I1_initial[p] = I1;
        particles.J2_initial[p] = J2;

        std::cout << "J2 = " << J2 << std::endl;

        // 3 cases of Return Mapping:

        // 1. perfect elasticity (f <= 0)

        T temp_criterion = f_dp(J2, p_I1, particles.k[p], eta);


        std::cout << "f Drucker-Prager" << std::endl;
        std::cout << temp_criterion << std::endl;

        int count = 0;

        if (temp_criterion <= 0.) {
	  // elastic case
	  particles.sigma[p] = (*stress_trial)[p];

        } else if (p_I1 > particles.k[p]/eta) {
	  /* projection to tip here */

	  particles.sigma[p] = TM2::Zero();

	  // particles.sigma[p] = 1./2.*( particles.k[p]/eta - particles.sigma_zz[p]) * TM2::Identity();

	  // T b = 0.5*((3.*particles.k[p]/eta - particles.sigma_zz[p]) + sqrt(-3.*particles.k[p]/eta*particles.k[p]/eta + 6. * particles.k[p]/eta * particles.sigma_zz[p] + particles.sigma_zz[p]*particles.sigma_zz[p] ));
	  // T b = 0.5*(3.*particles.k[p]/eta - particles.sigma_zz[p] - sqrt(3.)*abs(particles.k[p/eta - particles.sigma_zz[p]]));
	  // T a = 3.*particles.k[p]/eta - particles.sigma_zz[p] - b;

	  // particles.sigma[p](0,0) = a;
	  // particles.sigma[p](1,1) = b;
	  // particles.sigma[p](0,1) = 0.;
	  // particles.sigma[p](1,0) = 0.;

	  std::cout << "Project to tip" << std::endl;
	  // I1 = particles.sigma[p].trace() + particles.sigma_zz[p];
	  // p_I1 = I1/3.;
	  // s = particles.sigma[p] - p_I1*TM2::Identity();
	  // s_zz = particles.sigma_zz[p] - p_I1;
	  // J2 = s(0,1)*s(0,1) - s(0,0)*s(1,1) - s(0,0)*s_zz - s(1,1)*s_zz + 1.e-15;

	  // std::cout << " I1 " << I1 <<  "   " << " J2 " << J2 << "   " << " f_dp " << f_dp(J2, p_I1, particles.k[p], eta) << std::endl;
	  count = -1;
	  // particles.sigma[p] = I1/2.*TM2::Identity();

        } else if (temp_criterion > 0.) {
	  // Return Mapping to smooth part with Drucker-Prager line
	  std::cout << "RMA" << std::endl;

            // while (f_dp(J2, p_I1, particles.k[p], eta) > 0.) {
	  for (int i = 0; i < 20; i++) {
	    // T Delta_Lambda = (G/sqrt(J2) * s.cwiseProduct(strain_increment[p]).sum() + 3. * eta * K * strain_trace) / (G + 9. * eta * eta_g * K + xi * Hp * xi);
	    // T Delta_Lambda = f_dp(J2, p_I1, particles.k[p], eta)/(G + 4. * eta * eta_g * K + xi * Hp * xi);

	    // plastic multiplier
	    T Delta_Lambda = f_dp(J2, p_I1, particles.k[p], eta) / (G_3D + eta * eta_g * K_3D + xi * Hp * xi);

	    // updates
	    (*stress_trial)[p] -= Delta_Lambda * (G_3D/sqrt(J2) * s + K_3D * eta_g * TM2::Identity());
	    particles.k[p] += Delta_Lambda * xi * Hp * xi;

	    // recompute plastic quantities
	    I1 = (*stress_trial)[p].trace() + particles.sigma_zz[p];
	    p_I1 = I1/3.;
	    s = (*stress_trial)[p] - p_I1*TM2::Identity();
	    s_zz = particles.sigma_zz[p] - p_I1;
	    J2 = s(0,1)*s(0,1) - s(0,0)*s(1,1) - s(0,0)*s_zz - s(1,1)*s_zz + 1.e-15;
	    ++count;
	    std::cout << f_dp(J2, p_I1, particles.k[p], eta) << std::endl;
	    if (f_dp(J2, p_I1, particles.k[p], eta) < 1.e-05) {
	      break;
	    }
	  }

	  std::cout << count << std::endl;
	  // std::cout << f_dp(J2, p_I1, particles.k[p], eta) << std::endl;

	  particles.sigma[p] = (*stress_trial)[p];
        }

	// recompute everything to save and display after projection on the yield curve
	particles.I1[p] = particles.sigma[p].trace() + particles.sigma_zz[p];// p' = -1./3.*I1;
	p_I1 = particles.I1[p]/3.;
	s = particles.sigma[p] - p_I1*TM2::Identity();
	s_zz = particles.sigma_zz[p] - p_I1;
	particles.J2[p] = s(0,1)*s(0,1) - s(0,0)*s(1,1) - s(0,0)*s_zz - s(1,1)*s_zz + 1.e-15;// q = sqrt(3./2.*J2);
	particles.k_final[p] = particles.k[p];
	particles.f[p] = f_dp(particles.J2[p], p_I1, particles.k[p], eta);

	particles.count[p] = count;

	std::cout << "f Drucker-Prager final = " << particles.f[p] << std::endl;

      }
    }
  } else if (plasticity_model_int == 2) { // Modified Cam Clay

    auto f_mcc = [](T &q, T &p_I1, T &beta, T &p_0, T &M) { return q*q*(1 + 2*beta) + M*M*(p_I1 + beta*p_0)*(p_I1 - p_0); };

    auto hardening = [](T &K, T &xhi, T &epsilon) {return K*sinh(xhi*epsilon); };
    auto hardening_prime = [](T &K, T &xhi, T &epsilon) {return 0.; };
    auto Psi_q = [](T &q, T &G, T &beta, T &gamma) {return q/(1. + 6.*G*gamma*(1. + 2.*beta)); };
    auto Psi_p = [](T &p, T &K, T &epsilon_vp, T &epsilon) {return p - K*(epsilon - epsilon_vp); };
    
    particles.sigma.resize(Np);    std::fill(particles.sigma.begin(),    particles.sigma.end(),    TM2::Zero());

    unsigned int c_plastic = 0; 

#pragma omp parallel num_threads(n_threads)
    {
      std::vector<TM2> particles_sigma_local(Np, TM2::Zero());
      std::vector<T> particles_sigma_zz_local(Np, 0.);
      // std::vector<TM2> particles_sigma_local(Np, TM2::Zero());
      // std::vector<T> particles_sigma_zz_local(Np, 0.);
      // std::vector<T> particles_d_epsilon_v_p_local(Np, 0.);
      // std::vector<T> particles_epsilon_v_p_cumulative_local(Np, 0.);
      // std::vector<T> particles_p0_local(Np, 0.);
      // unsigned int c_plastic_local = 0;
      
      #pragma omp for
      for (int p=0; p<Np; p++) {
	particles_sigma_zz_local[p] = - 0.5*particles.rho[p]*particles.h[p]*gravity*cos(particles.theta[p]);

	// first invariant of Cauchy stress tensor
	T I1 = (*stress_trial)[p].trace();
	T p_I1 = -I1/2.;

	// deviatoric Cauchy stress tensor
	TM2 s = (*stress_trial)[p] + p_I1*TM2::Identity();

	// second invariant of deviatoric stress tensor
	// T J2 = (1/6) * ( (s(0,0) - s(1,1))*(s(0,0) - s(1,1)) + s(0,0)*s(0,0) +s(1,1)*s(1,1) + s(0,1)*s(0,1) ) + 1e-15;
	// T q = sqrt(3*J2);
	T J2 = s(0, 0)*s(0, 0) + s(1, 0)*s(1, 0) + s(0, 1)*s(0, 1) + s(1, 1)*s(1, 1);
	T q = sqrt((3./2.)*J2);
      
	particles.I1_initial[p] = I1;
	particles.J2_initial[p] = J2;

	T temp_criterion = f_mcc(q, p_I1, beta, particles.p0[p], M);

	if (temp_criterion <= 0. && particles.p0[p] > 1/2.) {
	  // pur elastic case
	  particles_sigma_local[p] = (*stress_trial)[p];
	  particles.I1[p] = I1;
	  particles.J2[p] = J2;

        } else  {
	  /* projection to tip here */
	  //std::cout << "Simplified RMA" << std::endl;
	  particles_sigma_local[p] = TM2::Zero();
	  particles.I1[p] = 0.;
	  particles.J2[p] = 0.;
	  particles.p0[p] = 0.;
    particles.count[p] = 0.;
	}


      // 	if (temp_criterion <= 0.) {
      // 	  // pur elastic case
      // 	  particles_sigma_local[p] = (*stress_trial)[p];
      // 	  particles.I1[p] = I1;
      // 	  particles.J2[p] = J2;

      //   } else  {
      // 	  c_plastic_local++;
      // 	  // std::cout << p_I1 << std::endl;
      // 	  // std::cout << q << std::endl;
      // 	  // std::cout << "RMA" << std::endl;
      // 	  // We choose initial for Delta gamma and epsilon plastic volumetric
      // 	  T gamma = 0.;
      // 	  T epsilon = 0.;
      // 	  TV2 diff;
      // 	  // We define jacobian and function
      // 	  TV2 Psi_MCC = {1., 1.};
      // 	  TM2 J_Psi;
      // 	  int c = 0;
      // 	  T epsilon_c = particles.epsilon_v_p_cumulative[p];

      //   while (Psi_MCC.norm() > threshold_MCC){
      // 	  T p_aux = Psi_p(p_I1, K, epsilon_c, epsilon);
      // 	  T q_square = Psi_q(q, G, beta, gamma)*Psi_q(q, G, beta, gamma);
      //     T p0 = particles.p0[p];
      //     // T p0 =  hardening(K, xi, epsilon_c);
      //     Psi_MCC(0) = (1. + 2.*beta)*q_square + M*M*(p_aux + beta*p0)*(p_aux - p0);

      //     Psi_MCC(1) = epsilon - epsilon_c - 2.*gamma*M*M*(p_aux + (beta-1)/2.*p0);

      //     J_Psi(0, 0) = -12.*G*(1.+2.*beta)*(1.+2.*beta)/(1. + 6.*G*gamma*(1.+2.*beta)) * q_square;
      //     J_Psi(0, 1) = M*M*((-K + beta*hardening_prime(K, xi, epsilon))*(p_aux - p0) + (p_aux + beta*p0)*(-K - hardening_prime(K, xi, epsilon)));

      //     J_Psi(1, 0) = -2*M*M*(p_aux + (beta-1)/2.*p0);

      //     J_Psi(1, 1) = 1. - 2*gamma*M*M*(-K + (beta-1.)/2.*hardening_prime(K, xi, epsilon));

      //     diff = - J_Psi.householderQr().solve(Psi_MCC);

      //     gamma +=diff(0);
      //     epsilon += diff(1);
      //     c++;
      //   }
      //   //std::cout << c << std::endl;
      //   p_I1 = Psi_p(p_I1, K, epsilon_c, epsilon);
      //   particles_d_epsilon_v_p_local[p] = epsilon - epsilon_c;
      //   particles_epsilon_v_p_cumulative_local[p] = epsilon;
      //   s = s/(1. + 6*G*gamma*(1. + 2.*beta));

      //   particles_sigma_local[p] = s + p_I1*TM2::Identity();
      //   T J2 = ( s(0,0)*s(0,0) + s(1,1)*s(1,1) + s(1, 0)*s(1, 0) + s(0,1)*s(0,1) ) + 1e-15;
      //   T q = sqrt(3*J2);

      //   particles.I1[p] = p_I1;
      //   particles.J2[p] = q;

      //   T final_criterium = f_mcc(q, p_I1, beta, particles.p0[p], M);
      //   //std::cout << "final criterium : ";
      //   //std::cout << final_criterium << std::endl;
      //   //std::cout << particles.p0[p] << std::endl;
      //   particles_p0_local[p] = std::max(hardening(K, xi, particles.epsilon_v_p_cumulative[p]), 0.);



      // 	}
       } // end for p

     // particles.d_epsilon_v_p.resize(Np); fill(particles.d_epsilon_v_p.begin(), particles.d_epsilon_v_p.end(), 0.);
     // //particles.epsilon_v_p_cumulative.resize(Np); fill(particles.epsilon_v_p_cumulative.begin(), particles.epsilon_v_p_cumulative.end(), 0.);
     //  particles.p0.resize(Np); fill( particles.p0.begin(), particles.p0.end(), 0. );

#pragma omp critical
      {
      	for (int p = 0; p < Np; p++) {
	  particles.sigma[p]    += particles_sigma_local[p];
	  particles.sigma_zz[p] += particles_sigma_zz_local[p];
	  // particles.sigma[p]    += particles_sigma_local[p];
          // particles.sigma_zz[p] += particles_sigma_zz_local[p];
          // particles.d_epsilon_v_p[p] += particles_d_epsilon_v_p_local[p];
          // particles.epsilon_v_p_cumulative[p] += particles_epsilon_v_p_cumulative_local[p];
          // particles.p0[p] += particles_p0_local[p];
	}
	// c_plastic += c_plastic_local;
      } // end omp critical
    }// end parallel

    // std::cout << "Number of Plastic Particles : ";
    // std::cout << c_plastic << std::endl;

  } else if (plasticity_model_int == 0) { // No Plasticity

    particles.sigma.resize(Np);    std::fill(particles.sigma.begin(),    particles.sigma.end(),    TM2::Zero());
    particles.sigma = (*stress_trial);
  }
};
