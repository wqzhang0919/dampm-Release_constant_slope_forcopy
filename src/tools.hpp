#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

/// FLOAT OR DOUBLE //// // DOUBLE !!!!!!!!!!
typedef double T;
////////////////////////
typedef Eigen::Matrix<T, 2, 2> TM2;
typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMX;
typedef Eigen::Matrix<T, 2, 1> TV2;
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVX;
typedef Eigen::Array<T, 2, 2> TA2;
////////////////////////

enum PlateType { top, bottom, left, right}; // need to define a Plate to be the slope?
enum ElasticModel { StvkWithHencky, NeoHookean };
enum PlasticModel { NoPlasticity, VonMises, DPSimpleSoft };
enum BoundaryCondition { STICKY, SLIP };

#define CUBICSPLINES
// #define QUADRATICSPLINES

///////////////////// TOOLS ////////////////////////

template <typename T>
void debug(T in){
  std::cout << in << std::endl;
}
template <typename T, typename U>
void debug(T in1, U in2){
  std::cout << in1 << in2 << std::endl;
}
template <typename T, typename U, typename V>
void debug(T in1, U in2, V in3){
  std::cout << in1 << in2 << in3 << std::endl;
}
template <typename T, typename U, typename V, typename W>
void debug(T in1, U in2, V in3, W in4){
  std::cout << in1 << in2 << in3 << in4 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X>
void debug(T in1, U in2, V in3, W in4, X in5){
  std::cout << in1 << in2 << in3 << in4 << in5 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X, typename Y>
void debug(T in1, U in2, V in3, W in4, X in5, Y in6){
  std::cout << in1 << in2 << in3 << in4 << in5 << in6 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X, typename Y, typename Z>
void debug(T in1, U in2, V in3, W in4, X in5, Y in6, Z in7){
  std::cout << in1 << in2 << in3 << in4 << in5 << in6 << in7 << std::endl;
}

inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

inline T selfDoubleDot(const TM2& A){
    T out = A(0,0)*A(0,0) + A(1,1)*A(1,1) + A(0,1)*A(0,1) + A(1,0)*A(1,0);
    return out;
}

void load_array(TVX& array_, unsigned int n_cols, std::string file_name);

std::vector<T> linspace(T a, T b, size_t N);

#ifdef CUBICSPLINES

    inline T N(const T &u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return 0.5 * uabs*uabs*uabs - uabs*uabs + 0.6666666666666666666666666666666666666666666666666;
        } else if (uabs < 2.0){
            return 0.1666666666666666666666666666666666666666666 * (2.0 - uabs) * (2.0 - uabs) * (2.0 - uabs);
        } else {
            return 0;
        }
    }

    inline T dNdu(const T &u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return u * (1.5 * uabs - 2.0);
        } else if (uabs < 2.0){
            return -0.5 * sgn(u) * (2.0 - uabs) * (2.0 - uabs);
        } else {
            return 0;
        }
    }

    inline T d2Ndu2(const T &u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return (3.0 * uabs - 2.0);
        } else if (uabs < 2.0){
            return (2.0 - uabs);
        } else {
            return 0;
        }
    }

#else // QUADRATICSPLINES // QUADRATIC SPLINES

    inline T N(const T &u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return 0.75 - uabs * uabs;
        } else if (uabs < 1.5){
            return 0.5 * (1.5 - uabs) * (1.5 - uabs);
        } else {
            return 0;
    	}
    }

    inline T dNdu(const T &u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return (-2*u);
        } else if (uabs < 1.5){
            return (u - 1.5*sgn(u));
        } else {
            return 0;
    	}
    }

    inline T d2Ndu2(const T &u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return -2.0;
        } else if (uabs < 1.5){
            return 1;
        } else {
            return 0;
    	}
    }

// #else //Linear splines

//     inline T N(const T &u){
//         T uabs = std::abs(u);
//         if (uabs < 1.){
//             return 1. - uabs;
//         } else {
//             return 0;
//     	}
//     }

//     inline T dNdu(const T &u){
//         T uabs = std::abs(u);
//         if (uabs < 1.){
//             return (T(0) < u) - (u > T(0));
//         } else {
//             return 0;
//     	}
//     }

//     inline T d2Ndu2(const T &u){
//             return 0;
//     }

#endif




inline T wip(const T &xp, const T &yp, const T &xi, const T &yi, const T &one_over_h){
    return N( (xp - xi) * one_over_h ) * N( (yp - yi) * one_over_h );
}

inline TV2 grad_wip(const T &xp, const T &yp, const T &xi, const T &yi, const T &one_over_h){
    TV2 out;
    out << dNdu((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) * one_over_h,
           dNdu((yp - yi) * one_over_h) * N((xp - xi) * one_over_h) * one_over_h;
    return out;
}

inline T laplace_wip(const T &xp, const T &yp, const T &xi, const T &yi, const T &one_over_h, const T &one_over_h_square){
    T term1 = d2Ndu2((xp - xi) * one_over_h) * N((yp - yi) * one_over_h);
    T term2 = d2Ndu2((yp - yi) * one_over_h) * N((xp - xi) * one_over_h);
    return ( term1 + term2 ) * one_over_h_square;
}

// grid numeration indices :
  inline unsigned int ind(const unsigned int & i, const unsigned int & j, const unsigned int & Ny){
      return i*Ny + j; // 2D
      // return (i*Ny + j) * Nz + k; // 3D
  }

// angles calculation
inline TV2 angles(const TV2 & part_terrain_gradients) { // computes slope angle and azimut of terrain gradients at a certain p particle
    TV2 out;
    T azimut;
    T slope = atan(part_terrain_gradients(0)); // quite strange init

    if (part_terrain_gradients(1) > 0.) {
        azimut = M_PI_2;
    } else {
        azimut = -M_PI_2;
    }

    if (part_terrain_gradients(0) != 0.) {
        azimut = atan(part_terrain_gradients(1)/part_terrain_gradients(0));
    }

    if (azimut != 0.) {
        slope = atan(part_terrain_gradients(1)/abs(sin(azimut)));
    }

    out << slope,
        azimut;

    return out;
}





#endif  // TOOLS_HPP
