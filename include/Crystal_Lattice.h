// Author: Cedric Schoonen, March 2025
// Contact: cedric.schoonen1@gmail.com

// Constructs a crystal lattice, orienting a specified plane orthogonally to the z axis
// The scale is such that interatomic distance between nearest neighbors is 1

#ifndef __SCHOONEN_CRYSTAL_LATTICE_
#define __SCHOONEN_CRYSTAL_LATTICE_

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>

using namespace std;


class Crystal_Lattice
{
  public:
  Crystal_Lattice(string lattice_name, string lattice_orientation,
                   int num_cells_x=1, int num_cells_y=1, int num_cells_z=1)
  {
    lattice_name_ = lattice_name;
    lattice_orientation_ = lattice_orientation;
    
    initialize_unit_cell();
    
    num_cells_[0] = num_cells_x;
    num_cells_[1] = num_cells_y;
    num_cells_[2] = num_cells_z;
    
    if (num_cells_x<=0 || num_cells_y<=0 || num_cells_z<=0)
      throw runtime_error("Error in Crystal_Lattice.h: The number of copies of the unit cell must be a strictly positive integer!");
    
    make_copies_of_the_unit_cell();
  }
  
  ~Crystal_Lattice() {}

  int get_num_cells(int i) const {return num_cells_[i];}
  double get_lattice_length(int i) const {return L_[i];}
  
  // Equal proportions in all directions, interatomic distance of 1 (between nearest neighbors)
  vector<vector<double>> get_atoms() const { return atoms_;}
  
  // Equal proportions in all directions, interatomic distances scaled from 1 to dnn
  vector<vector<double>> get_atoms_scaled(double dnn) const;
  
  // Lattice streched to match the given dimensions, inconsistent distances between original nearest neighbors
  vector<vector<double>> get_atoms_scaled(double Lx, double Ly, double Lz) const;
  
  void export_to_xyz_file(string filename="crystal_lattice.xyz") const;

  protected:
  void initialize_unit_cell();
  void make_copies_of_the_unit_cell();
  
  vector<vector<double>> atoms_; // stores lattice site positions
  
  string lattice_name_;        // accepts either BCC, FCC or HCP
  string lattice_orientation_; // accepts either 001, 110 or 111 for BCC/FCC, or 001,010,100 for HCP 
  
  int num_cells_[3] = {1,1,1}; // number of copies of the default rectangular unit cell
  double L_[3] = {0,0,0};      // side lengths of the lattice (assuming dnn=1)
};


#endif
