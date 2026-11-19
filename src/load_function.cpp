#include "simulation.hpp"
#include "tools.hpp"
using namespace std;

void Simulation::read_height(string extra){

  std::cout << "Begin Reading Height" << std::endl;
  // we first create un h vector  then we will rearrange it to match particles.h
  std::vector<T> h_vector;
  h_vector.resize(Np); std::fill(h_vector.begin(), h_vector.end(), 0.);
  std::string Dir = "/home/guillet/Work/Francis/dampm/src/";
  unsigned int index = 0;

  fstream Myfile;
  Myfile.open(Dir + height_file, ios::in);
  if (!Myfile) {
		std::cout << "Height file " << height_file << " not found, exit code" << std::endl;
    abort();
	} else {
		
    T h;
		while (1) {
			Myfile >> h;
      // std::cout << h << std::endl;
      if (index > Np-1) {std::cout << "Wrong size height file" << std::endl; abort();}
      h_vector[index] = h;
      index++;
			if (Myfile.eof()) {break;}		
		}
  // std::cout << "end truc" << std::endl;
	}
	Myfile.close();

  //at the end of the loop index must be equal at Np
  if (index < Np) {std::cout << "Wrong size height file" << std::endl; abort();}
  //Now we have to link h_vector to particles.h
  std::vector<T> col_i;
  unsigned int p_init;
  std::cout << "Ny : " << Ny << std::endl;

  // UNOPTIMIZED BUT MEAN TO CHANGED IN THE FUTURE
  for(int i=0; i<Nx; i++){
    col_i.resize(4*Ny); std::fill(col_i.begin(), col_i.end(), 0.);
    for(int j=0;j<Ny; j++){
      col_i[4*(Ny-j-1) + 3] = h_vector[4*j*Nx + 2*i +1];
      col_i[4*(Ny-j-1) + 2] = h_vector[4*j*Nx + 2*i];
      col_i[4*(Ny-j-1) + 1] = h_vector[4*j*Nx + 2*Nx + 2*i +1];
      col_i[4*(Ny-j-1)]     = h_vector[4*j*Nx + 2*Nx + 2*i];
    }
    p_init = i*4*Ny;
    for(int j=0; j<Ny; j++){
      particles.h[p_init + 4*j]     = col_i[4*j];
      particles.h[p_init + 4*j + 1] = col_i[4*j+3];
      particles.h[p_init + 4*j + 2] = col_i[4*j+2];
      particles.h[p_init + 4*j + 3] = col_i[4*j+1];
    }
  }

  std::cout <<"End Reading Height" << std::endl;

}

