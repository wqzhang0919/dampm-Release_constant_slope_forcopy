// #include "simulation.hpp"
#include "tools.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <math.h>
#include <vector>

using namespace std;



class Poisson_disk{
    private:
    int k = 30; //Number maximum of try to find a new particles
    T r; //Minimum distance between particles
    T Lx; //Size of the surface
    T Ly; //Size of the surface
    T a; //Cell side length
    int Nx; //Number of cells in the x direction
    int Ny; //Number of cells in the y direction
    vector<TV2> samples; //Result
    vector<int> cells; //Table to save the index of the particles in cells, equal to -1 if empty
    


    public:
    Poisson_disk(T dx, T L1, T L2, int ppc){
        r = sqrt(15./24.*dx*dx/ppc); //Minimum distance between particles to respect the numbers of particles per cell, 15/24 is an empirical coefficient (better than 2/3 of Jixie)
        Lx = L1; 
        Ly = L2;
        a = r/sqrt(2); //To be sure that it is not possible to have more than one particle in a cell
        Nx = int(Lx/a)+1;
        Ny = int(Ly/a)+1;

        cells.resize(Nx*Ny);    
        for(int i = 0; i < Nx; i++){
            for(int j = 0; j < Ny; j++){
                cells[i*Ny+j] = -1; //Initialization at -1 for all the cells means that there is no particle
            }
        }
    }


    private:
    vector<int> get_cell_coords(TV2 pt){
    /*Get the coordinates of the cell that pt = (x,y) falls in.*/

        vector<int> res;
        res.resize(2);

        res[0] = floor(pt(0)/a);
        res[1] = floor(pt(1)/a);

        return res;
    }


    private:
    vector<int> get_neighbours(vector<int> coords){
    /*Return the indexes of points in cells neighbouring cell at coords.

    For the cell at coords = (x,y), return the indexes of points in the cells
    with neighbouring coordinates illustrated below: ie those cells that could 
    contain points closer than r.
                                     ooo
                                    ooooo
                                    ooXoo
                                    ooooo
                                     ooo
    */



        int dxdy[21][2] = {{-1,-2},{0,-2},{1,-2},{-2,-1},{-1,-1},{0,-1},{1,-1},{2,-1},{-2,0},{-1,0},{1,0},{2,0},{-2,1},{-1,1},{0,1},{1,1},{2,1},{-1,2},{0,2},{1,2},{0,0}};
        vector<int> neighbours;

        for(int i = 0; i<21; i++){ //Test for all the neighbours

            vector<int> neighbour_coords;
            neighbour_coords.resize(2);

            neighbour_coords[0] = coords[0]+dxdy[i][0]; neighbour_coords[1] = coords[1]+dxdy[i][1];
            if (neighbour_coords[0]>=0 && neighbour_coords[0] <Nx && neighbour_coords[1]>=0 && neighbour_coords[1]<Ny) { //Test if we are not off the grid
                int neighbour_cell = cells[neighbour_coords[0]*Ny+neighbour_coords[1]]; 

                if(neighbour_cell != -1){//Test if there is a particle in the neighbour cell
                    neighbours.push_back(neighbour_cell);
                }

            }
        }
        return neighbours;
    }

    private:
    bool point_valid(TV2 pt){
    /*Is pt a valid point to emit as a sample?

    It must be no closer than r from any other point: check the cells in its
    immediate neighbourhood.*/

        vector<int> cell_coords = get_cell_coords(pt);
        vector<int> neighbours = get_neighbours(cell_coords);

        for(int i=0; i<neighbours.size();i++){//Check for all the neighbours
            TV2 nearby_pt = samples[neighbours[i]];
            T distance2 = (nearby_pt[0]-pt[0])*(nearby_pt[0]-pt[0])+(nearby_pt[1]-pt[1])*(nearby_pt[1]-pt[1]);

            if(distance2 <r*r){
                return false;
            }
        }
        return true;
    }


    private:
    TV2 get_point(TV2 refpt){
    /*Try to find a candidate point relative to refpt to emit in the sample.

    We draw up to k points from the annulus of inner radius r, outer radius 2r
    around the reference point, refpt. If none of them are suitable (because
    they're too close to existing points in the sample), return False.
    Otherwise, return the pt.*/

        for(int i=0;i<k;i++){
            T rho = float(std::rand())/float(RAND_MAX)*r + r;
            T theta = float(std::rand())/float(RAND_MAX)*2*M_PI;
            TV2 pt = {refpt[0] + rho*cos(theta), refpt[1] + rho*sin(theta)};
            if(pt[0]>=0 && pt[0]<= Lx && pt[1]>=0 && pt[1]<=Ly){
                if(point_valid(pt)){
                    return pt;
                }
            }
        }
        return {-1,-1}; //If no point found
    }

    public:
    vector<TV2> poisson_disk_sampling(){
        /*Return the list of the particles position based on the algorithm Fast Poisson Disk Sampling in Arbitrary Dimensions by Robert Bridson*/

        TV2 pt = {float(std::rand())/float(RAND_MAX)*Lx,float(std::rand())/float(RAND_MAX)*Ly}; //Pick a random point to start with
        samples.push_back(pt);
        vector<int> cell_coords = get_cell_coords(pt);
        
        cells[cell_coords[0]*Ny+cell_coords[1]] = 0; //Our first sample is indexed at 0 in the samples list...
        vector<int> active = {0}; //... and it is active, in the sense that we're going to look for more points in its neighbourhood.

        int nsamples = 0; //The index of the current sample
        while (!active.empty()) { //As long as there are points in the active list, keep trying to find samples.

            int idx = trunc(float(std::rand())/float(RAND_MAX)*active.size()); //Choose a random "reference" point from the active list.
            TV2 refpt = samples[active[idx]];

            pt = get_point(refpt); //Try to pick a new point relative to the reference point.

            if(pt[0] == -1 && pt[1] == -1){
                active.erase(active.begin()+idx); //We had to give up looking for valid points near refpt, so remove it from the list of "active" points.
            }
            else {
                //Point pt is valid: add it to the samples list and mark it as active
                samples.push_back(pt);
                nsamples+=1;

                active.push_back(nsamples);

                cell_coords = get_cell_coords(pt);
                cells[cell_coords[0]*Ny+cell_coords[1]] = nsamples;
            }
        }
        return samples;
    }
};