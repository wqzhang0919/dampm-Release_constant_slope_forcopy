#ifndef OBJECT_HPP
#define OBJECT_HPP

#include "tools.hpp"
#include <string>

class InfinitePlate{
public:
  InfinitePlate() : x_object(0.0), y_object(0.0), vx_object(0.0), vy_object(0.0), plate_type(bottom), name("noname") {}
  InfinitePlate(T x_object, T y_object, T vx_object, T vy_object, PlateType plate_type, std::string name) : x_object(x_object), y_object(y_object), vx_object(vx_object), vy_object(vy_object), plate_type(plate_type), name(name) {}

  bool inside(T x, T y){
      return distance(x, y) <= 0; // inside if dist is negative
  }

// add slope plate here with angle theta! 
  T distance(T x, T y){
      if (plate_type == bottom)
        return (y - y_object);
      else if (plate_type == top)
          return (y_object - y);
      else if (plate_type == left)
          return (x - x_object);
      else // right
          return (x_object - x);
  }

  void move(T dt){
      x_object += dt * vx_object;
      y_object += dt * vy_object;
  }

  T x_object;
  T y_object;

  T vx_object;
  T vy_object;

  PlateType plate_type;
  std::string name;

};




#endif  // OBJECT_HPP
