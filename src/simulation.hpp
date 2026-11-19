#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>

#include <vector>
#include <chrono>
#include "tools.hpp"
#include "object.hpp"
#include "data_arrays.hpp"
#include "timer.hpp"
#include "poisson_disk_sampling.cpp"

// #define OMP //comment out if not using omp

class Simulation{
public:
  //Simulation() : current_time_step(0), exit(0) {};
  Simulation();
  void initialize(T density);
  void initialize_particles(unsigned int const& NloopX, unsigned int const& NloopY, unsigned int const& ppc, T const& material_density, T const& slope_angle);
  void simulate();
  void saveParticle(std::string extra = "");
  void saveParticleData(std::string extra = "");
  void saveParticleData_light(std::string extra = "");
  void saveParticleData_PST(std::string extra = "");
  void saveParticleData_2D(std::string extra = "");
  void saveGridData(std::string extra = "");
  void createDirectory();
  void read_config(std::string extra = "");
  void read_height(std::string extra = "");
  TV2 terrain_gradients(const int & p);

  bool height_loading;
  bool config_loading;

  std::string height_file;
  std::string config_file;
  std::string sim_name;
  std::string directory;
  std::string sim_type;
  std::string plasticity_model;
  std::string saving_type;
  int plasticity_model_int;
  int sim_type_int;
  bool compressibility;
  bool parallelised;
  bool poisson_disk;

  unsigned int dim; // always 2 here
  unsigned int n_threads; // to check and remove
  unsigned int current_time_step;
  unsigned int frame;
  unsigned int end_frame;
  int exit;
  T time;
  T final_time;
  T frame_dt;
  T dt;
  T dt_max;
  T loading_time;
  T safety_factor;
  T velocity_PST;

  T flip_ratio;
  T wave_speed;
  T cfl;
  T dx;
  // T L;
  T L1;
  T L2;
  T rho;
  T rho_water;
  T void_ratio;
  // TV2 gravity;
  T gravity;

  T amplitude;

      T fluidization;
    T basal_frictional_angle;

  // Particle data
  unsigned int Np;
  // T particle_mass;
  // T particle_area; // initial particle volume V0
  Particles particles;

  // Grid data
  unsigned int Nx, Ny, Nz;
  Grid grid;
  // inline unsigned int ind(const unsigned int & i, const unsigned int & j, const unsigned int & Ny){
  //     return i*Ny + j; // 2D
  //     // return (i*Ny + j) * Nz + k; // 3D
  // }

  // Elastoplasticity
  T mu;
  T lambda;
  // T K;
  T phi; // sidewall frictional angle
  T psi; // dilatancy angle
  T cohesion;
  T E;
  T nu;
  T threshold_MCC;


  //Weak Layer
  T p0_Initial;
  T Initial_cohesion_WL;
  T Thickness_WL;
  T Stiffness_WL;
  T Shear_modulus_WL;
  T Shear_strength;

  // Objects
  T friction;
  BoundaryCondition boundary_condition;
  std::vector<InfinitePlate> objects;

  // Runtime measurements
  T runtime_p2g;
  T runtime_g2p;
  T runtime_g2p_1;
  T runtime_g2p_2;
  T runtime_g2p_3;
  T runtime_comp_for;
  T runtime_elasticity;
  T runtime_plasticity;

  // Elasticity quantities, expressed in 2D (plane stress assumption)
  T G; // shear modulus
  T K; // bulk modulus

  // Elasticity quantities, expressed in 3D (for plasticity with sigma_zz force)
  T G_3D;
  T K_3D;

  // Plasticity quantities (Drucker Prager)
  T eta; // slope of Drucker Prager yield criterion
  T eta_g; // slope of the plastic potential g
  T xi; // parameter of hardening k
  T Hp; // plastic hardening modulus;

  // Plasticity quantities (MCC)
  T beta;
  T p_0;
  T M;

  // Precomputations
  T one_over_dx;
  T one_over_dx_square;

  // Functions
  void advanceStep();
  void updateDt();
  void remesh();

  void P2G();
  void P2G_dampm();
  void P2G_dampm_omp();

  void G2P();
  void G2P_dampm();
  void G2P_dampm_omp();
  void G2P_dampm_omp_1();
  void G2P_dampm_omp_2();
  void G2P_dampm_omp_3();

  void positionUpdate();

  void boundaryCollision(T xi, T yi, TV2& vi);
  // void boundaryCorrection(T xi, T yi, T& vxi, T& vyi);
//   void overwriteGridVelocity(T xi, T yi, TV2& vi);
//   void calculateMomentumOnParticles();
//   void calculateMomentumOnGrid();
//   void calculateMassConservation();

  void addExternalParticleGravity();
  std::pair<TMX, TMX> createExternalGridGravity();

  void compute_forces();
  void compute_forces_flow();
  void compute_forces_flow_omp();
  void compute_forces_release();
  void compute_forces_release_loading();
  void compute_forces_release_loading_omp();
  void compute_forces_release_omp();

  void compute_accel();
  void update();
  void update_position();
  // void strain_stress();
  std::unique_ptr< std::vector<TM2>> elasticity();
  std::unique_ptr< std::vector<TM2>> elasticity_omp();
  void plasticity(std::unique_ptr< std::vector<TM2>> & stress_trial);
  void plasticity_omp(std::unique_ptr< std::vector<TM2>> & stress_trial);

  //std::vector<TM2> elasticity();
  //std::vector<TM2> elasticity_omp();
  //void plasticity(std::vector<TM2> stress_trial);
  //void plasticity_omp(std::vector<TM2> stress_trial);

  void update_cohesion();
  void update_angle();

  void update_height();
}; // END Simulation class


#endif  // SIMULATION_HPP
