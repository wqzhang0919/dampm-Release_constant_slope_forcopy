#ifndef DATAARRAYS_HPP
#define DATAARRAYS_HPP

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "tools.hpp"

class Particles{
public:
  Particles(unsigned int Np = 400){
      x.resize(Np);    std::fill( x.begin(),    x.end(),    TV2::Zero() ); // position of particles in 2D
      v.resize(Np);    std::fill( v.begin(),    v.end(),    TV2::Zero() ); // velocity of particles in 2D
      a.resize(Np);    std::fill( a.begin(),    a.end(),    TV2::Zero() );
      
      h.resize(Np);   std::fill( h.begin(),   h.end(),   0.0); // height of particles 
      A.resize(Np);   std::fill( A.begin(),   A.end(),   0.0);
      x0.resize(Np); std::fill( x0.begin(), x0.end(), TV2::Zero() ); //Initial positions
      theta.resize(Np); std::fill( theta.begin(), theta.end(), 0.);//35./180.*M_PI); //angle for each particule (constant)
      rho.resize(Np); std::fill(rho.begin(), rho.end(), 0.);
      surfacic_mass.resize(Np); std::fill(surfacic_mass.begin(), surfacic_mass.end(), 0.);
      I1.resize(Np); std::fill(I1.begin(),  I1.end(), 0.);
      J2.resize(Np); std::fill(J2.begin(),  J2.end(), 0.);
      I1_initial.resize(Np); std::fill(I1_initial.begin(),  I1_initial.end(), 0.);
      J2_initial.resize(Np); std::fill(J2_initial.begin(),  J2_initial.end(), 0.);
      d_epsilon_v_p.resize(Np); std::fill(d_epsilon_v_p.begin(), d_epsilon_v_p.end(), 0.);
      epsilon_v_p_cumulative.resize(Np); std::fill(epsilon_v_p_cumulative.begin(), epsilon_v_p_cumulative.end(), 0.);
      k.resize(Np); std::fill(k.begin(),    k.end(), 0.0);
      k_initial.resize(Np); std::fill(k_initial.begin(),    k_initial.end(), 0.);
      k_final.resize(Np); std::fill(k_final.begin(),        k_final.end(),   0.);
      f.resize(Np); std::fill(f.begin(),        f.end(), 0.);
      old_position.resize(Np); std::fill(old_position.begin(),  old_position.end(), TV2::Zero());
      mass.resize(Np); std::fill(mass.begin(),  mass.end(), 0.0);

      sigma.resize(Np); std::fill(sigma.begin(), sigma.end(), TM2::Identity());

      sigma_xz.resize(Np); std::fill(sigma_xz.begin(),  sigma_xz.end(), 0.0);
      sigma_yz.resize(Np); std::fill(sigma_yz.begin(),  sigma_yz.end(), 0.0);
      sigma_zz.resize(Np); std::fill(sigma_zz.begin(),  sigma_zz.end(), 0.0);

      pressure.resize(Np); std::fill(pressure.begin(),  pressure.end(), 0.0);

      count.resize(Np); std::fill(count.begin(),        count.end(), 0.);

    //For the release, we keep in memory
      cohesion_WL.resize(Np); std::fill(cohesion_WL.begin(), cohesion_WL.end(), 0.);
      pic.resize(Np); std::fill(pic.begin(), pic.end(), TV2::Zero() );
      flip.resize(Np); std::fill( flip.begin(), flip.end(),  TV2::Zero() );
      p0.resize(Np); std::fill( p0.begin(), p0.end(), 0. );

    // For APIC
      B.resize(Np); std::fill( B.begin(), B.end(), TM2::Zero() );

  }

  std::vector<TV2> x;
  std::vector<TV2> x0;
  std::vector<TV2> v;
  std::vector<TV2> a;

  std::vector<TV2> old_position;
  std::vector<T> h;
  std::vector<T> theta;
  std::vector<T> rho;
  std::vector<T> surfacic_mass;
  std::vector<T> I1;
  std::vector<T> J2;
  std::vector<T> I1_initial;
  std::vector<T> J2_initial;
  std::vector<T> k;
  std::vector<T> k_initial;
  std::vector<T> k_final;
  std::vector<T> f;
  std::vector<T> A;
  std::vector<T> mass;
  std::vector<T> d_epsilon_v_p;
  std::vector<T> epsilon_v_p_cumulative;
  std::vector<TM2> sigma;

  std::vector<T> sigma_xz;
  std::vector<T> sigma_yz;
  std::vector<T> sigma_zz;

  std::vector<T> pressure;

  std::vector<T> count;

  std::vector<TV2> flip;
  std::vector<TV2> pic;
  std::vector<T> cohesion_WL;
  std::vector<T> p0;

  // Apic quantities
  std::vector<TM2> B;
};

class Grid{
public:
    Grid(unsigned int const& Nx = 21, unsigned int const& Ny = 21){
        x.resize(Nx); std::fill( x.begin(), x.end(), 0.0 ); // x component of position of nodes
        y.resize(Ny); std::fill( y.begin(), y.end(), 0.0 ); // y component of position of nodes

        v.resize(Nx*Ny);    std::fill( v.begin(),    v.end(),    TV2::Zero() ); // velocity in 2D of nodes
        a.resize(Nx*Ny);    std::fill( a.begin(),    a.end(),    TV2::Zero() ); // accel in 2D 

        mass.resize(Nx*Ny); std::fill( mass.begin(), mass.end(), 0.0 ); // nodal masses
        F.resize(Nx*Ny);    std::fill(F.begin(),    F.end(),    TV2::Zero());

        z.resize(Nx*Ny);    std::fill( z.begin(),    z.end(),    0.);

    }
    std::vector<T> x;
    std::vector<T> y;
    std::vector<TV2> v;
    std::vector<TV2> a;
    std::vector<T> mass;
    std::vector<TV2> F;
    std::vector<T> z; // real elevation = z_0 + h    
    T xc; // corner (x component) // to delete
    T yc; // corner (y component)
};


#endif  // DATAARRAYS_HPP
