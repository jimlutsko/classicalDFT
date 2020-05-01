#ifndef __LUTSKO__POINT__
#define __LUTSKO__POINT__

/**
  *  @brief Point Class: holds x,y,z coordinates
  *
  *    This class is just a triple with members named X,Y and Z
  */  


class Point
{
 public:
  /**
  *   @brief  Default  constructor for Point: all corrdinates are set to zero
  *  
  *   @return nothing 
  */  
 Point() : x_(0.0),y_(0.0),z_(0.0) {};

  /**
  *   @brief  Constructor for Point: 
  *  
  *   @param  x is a coordinate
  *   @param  y is a coordinate
  *   @param  z is a coordinate
  *   @return nothing 
  */  
 Point(double x, double y, double z) : x_(x),y_(y),z_(z) {};

  /**
  *   @brief  Accessor for x coordinate
  *  
  *   @return x_ 
  */  
  double X() const {return x_;}

  /**
  *   @brief  Accessor for y coordinate
  *  
  *   @return y_
  */  
  double Y() const {return y_;}

  /**
  *   @brief  Accessor for z coordinate
  *  
  *   @return z_ 
  */  
  double Z() const {return z_;}

 protected:
  double x_;  ///< x coordinate
  double y_;  ///< y coordinate
  double z_;  ///< z coordinate
};

#endif // __LUTSKO__POINT__
