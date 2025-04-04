#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>

using namespace std;

#include "Crystal_Lattice.h"


// Here we assume side lengths of 1 for the unit cell
void scale_atoms(vector<vector<double>> &atoms, double L[3])
{
  for (int i=0; i<atoms.size(); i++)
  {
    atoms[i][0] *= L[0];
    atoms[i][1] *= L[1];
    atoms[i][2] *= L[2];
  }
}


// Here we assume side lengths of 1 for the unit cell
void init_atoms_regular_bcc(vector<vector<double>> &atoms)
{
  atoms = vector<vector<double>>(2, vector<double>(3,0.0));
  
  atoms[0][0] = 0.0;
  atoms[0][1] = 0.0;
  atoms[0][2] = 0.0;

  atoms[1][0] = 1.0/2;
  atoms[1][1] = 1.0/2;
  atoms[1][2] = 1.0/2;
}


// Here we assume side lengths of 1 for the unit cell
void init_atoms_regular_fcc(vector<vector<double>> &atoms)
{
  atoms = vector<vector<double>>(4, vector<double>(3,0.0));
  
  atoms[0][0] = 0.0;
  atoms[0][1] = 0.0;
  atoms[0][2] = 0.0;

  atoms[1][0] = 1.0/2;
  atoms[1][1] = 1.0/2;
  atoms[1][2] = 0.0;

  atoms[2][0] = 1.0/2;
  atoms[2][1] = 0.0;
  atoms[2][2] = 1.0/2;

  atoms[3][0] = 0.0;
  atoms[3][1] = 1.0/2;
  atoms[3][2] = 1.0/2;
}


// Here we assume side lengths of 1 for the unit cell
void init_atoms_hcp_001(vector<vector<double>> &atoms)
{
  atoms = vector<vector<double>>(4, vector<double>(3,0.0));

  atoms[0][0] = 0.0;
  atoms[0][1] = 0.0;
  atoms[0][2] = 0.0;

  atoms[1][0] = 1.0/2;
  atoms[1][1] = 1.0/2;
  atoms[1][2] = 0.0;

  atoms[2][0] = 0.0;
  atoms[2][1] = 2.0/6;
  atoms[2][2] = 1.0/2;

  atoms[3][0] = 1.0/2;
  atoms[3][1] = 5.0/6;
  atoms[3][2] = 1.0/2;
}


// Defines atomic positions in the unit cell, as well as the side lengths
// Here the unit of distance is the distance between nearest neighbors
void Crystal_Lattice::initialize_unit_cell()
{
  if ( lattice_name_=="BCC" && lattice_orientation_=="001" ||
       lattice_name_=="BCC" && lattice_orientation_=="010" ||
       lattice_name_=="BCC" && lattice_orientation_=="100" )
  {
    init_atoms_regular_bcc(atoms_);
    
    for (int i=0; i<3; i++) L_[i] = 2.0/sqrt(3);
    
    scale_atoms(atoms_, L_);
  }
  else if ( lattice_name_=="BCC" && lattice_orientation_=="110" ||
            lattice_name_=="BCC" && lattice_orientation_=="101" ||
            lattice_name_=="BCC" && lattice_orientation_=="011" )
  {
    init_atoms_regular_fcc(atoms_); // yes, same relative positions in the unit cell
    
    for (int i=0; i<3; i++) L_[i] = 2.0/sqrt(3);
    if (lattice_orientation_=="110") {L_[0] *= sqrt(2); L_[1] *= sqrt(2);}
    if (lattice_orientation_=="101") {L_[0] *= sqrt(2); L_[2] *= sqrt(2);}
    if (lattice_orientation_=="011") {L_[1] *= sqrt(2); L_[2] *= sqrt(2);}
    
    scale_atoms(atoms_, L_);
  }
  else if ( lattice_name_=="BCC" && lattice_orientation_=="111" )
  {
    atoms_ = vector<vector<double>>(12, vector<double>(3,0.0));
    
    atoms_[0][0] = 1.0/6;
    atoms_[0][1] = 1.0/2;
    atoms_[0][2] = 0.0;

    atoms_[1][0] = 0.0;
    atoms_[1][1] = 0.0;
    atoms_[1][2] = 1.0/3;

    atoms_[2][0] = 1.0/6;
    atoms_[2][1] = 1.0/2;
    atoms_[2][2] = 1.0/2;

    atoms_[3][0] = 0.0;
    atoms_[3][1] = 0.0;
    atoms_[3][2] = 5.0/6;

    atoms_[4][0] = 3.0/6;
    atoms_[4][1] = 1.0/2;
    atoms_[4][2] = 1.0/3;

    atoms_[5][0] = 1.0/3;
    atoms_[5][1] = 0.0;
    atoms_[5][2] = 2.0/3;

    atoms_[6][0] = 3.0/6;
    atoms_[6][1] = 1.0/2;
    atoms_[6][2] = 5.0/6;

    atoms_[7][0] = 1.0/3;
    atoms_[7][1] = 0.0;
    atoms_[7][2] = 1.0/6;

    atoms_[8][0] = 5.0/6;
    atoms_[8][1] = 1.0/2;
    atoms_[8][2] = 2.0/3;

    atoms_[9][0] = 2.0/3;
    atoms_[9][1] = 0.0;
    atoms_[9][2] = 0.0;

    atoms_[10][0] = 5.0/6;
    atoms_[10][1] = 1.0/2;
    atoms_[10][2] = 1.0/6;

    atoms_[11][0] = 2.0/3;
    atoms_[11][1] = 0.0;
    atoms_[11][2] = 1.0/2;
    
    L_[0] = 2.0*sqrt(2);
    L_[1] = 2.0*sqrt(2.0/3);
    L_[2] = 2.0;
    
    scale_atoms(atoms_, L_);
  }
  else if ( lattice_name_=="FCC" && lattice_orientation_=="001" ||
            lattice_name_=="FCC" && lattice_orientation_=="010" ||
            lattice_name_=="FCC" && lattice_orientation_=="100" )
  {
    init_atoms_regular_fcc(atoms_);
    
    for (int i=0; i<3; i++) L_[i] = sqrt(2);
    
    scale_atoms(atoms_, L_);
  }
  else if ( lattice_name_=="FCC" && lattice_orientation_=="110" ||
            lattice_name_=="FCC" && lattice_orientation_=="101" ||
            lattice_name_=="FCC" && lattice_orientation_=="011" )
  {
    init_atoms_regular_bcc(atoms_); // yes, same relative positions in the unit cell
    
    for (int i=0; i<3; i++) L_[i] = 1.0;
    if (lattice_orientation_=="110") {L_[2] *= sqrt(2);}
    if (lattice_orientation_=="101") {L_[1] *= sqrt(2);}
    if (lattice_orientation_=="011") {L_[0] *= sqrt(2);}
    
    scale_atoms(atoms_, L_);
  }
  else if ( lattice_name_=="FCC" && lattice_orientation_=="111" )
  {
    atoms_ = vector<vector<double>>(6, vector<double>(3,0.0));

    atoms_[0][0] = 0.0;
    atoms_[0][1] = 0.0;
    atoms_[0][2] = 0.0; //A

    atoms_[1][0] = 1.0/2;
    atoms_[1][1] = 1.0/2;
    atoms_[1][2] = 0.0; //A

    atoms_[2][0] = 1.0/2;
    atoms_[2][1] = 1.0/6;
    atoms_[2][2] = 1.0/3; //B

    atoms_[3][0] = 0.0;
    atoms_[3][1] = 4.0/6;
    atoms_[3][2] = 1.0/3; //B

    atoms_[4][0] = 0.0;
    atoms_[4][1] = 2.0/6;
    atoms_[4][2] = 2.0/3; //C

    atoms_[5][0] = 1.0/2;
    atoms_[5][1] = 5.0/6;
    atoms_[5][2] = 2.0/3; //C
    
    L_[0] = 1.0;
    L_[1] = sqrt(3);
    L_[2] = sqrt(6);
    
    scale_atoms(atoms_, L_);
  }
  else if ( lattice_name_=="HCP" && lattice_orientation_=="001" )
  {
    init_atoms_hcp_001(atoms_);
    
    L_[0] = 1.0;
    L_[1] = sqrt(3);
    L_[2] = sqrt(8.0/3);
    
    scale_atoms(atoms_, L_);
  }
  else if ( lattice_name_=="HCP" && lattice_orientation_=="010" )
  {
    init_atoms_hcp_001(atoms_);
    
    // rotate along x axis so that the old y coordinate becomes the new z
    for (int i=0; i<atoms_.size(); i++)
    {
      double temp = atoms_[i][2];
      atoms_[i][2] = atoms_[i][1];
      atoms_[i][1] = -temp;
    }
    
    L_[0] = 1.0;
    L_[1] = sqrt(8.0/3);
    L_[2] = sqrt(3);
    
    scale_atoms(atoms_, L_);
  }
  else if ( lattice_name_=="HCP" && lattice_orientation_=="100" )
  {
    init_atoms_hcp_001(atoms_);
    
    // rotate along x axis so that the old x coordinate becomes the new z
    for (int i=0; i<atoms_.size(); i++)
    {
      double temp = atoms_[i][2];
      atoms_[i][2] = atoms_[i][0];
      atoms_[i][0] = -temp;
    }
    
    L_[0] = sqrt(8.0/3);
    L_[1] = sqrt(3);
    L_[2] = 1.0;
    
    scale_atoms(atoms_, L_);
  }
  else
  {
    throw runtime_error("Error in Crystal_Lattice.h: Unknown combination of lattice name + orientation!");
  }
}


void Crystal_Lattice::make_copies_of_the_unit_cell()
{
  // Record a copy of the unit cell
  vector<vector<double>> atoms_unit_cell = atoms_;
  double L_unit_cell[3] = {L_[0],L_[1],L_[2]};
  
  // Scale lattice side lengths
  for (int i=0; i<3; i++) L_[i] *= num_cells_[i];
  
  // Copy atoms
  atoms_.clear();
  for (int ix=0; ix<num_cells_[0]; ix++)
  for (int iy=0; iy<num_cells_[1]; iy++)
  for (int iz=0; iz<num_cells_[2]; iz++)
  for (int a=0; a<atoms_unit_cell.size(); a++)
  {
    vector<double> new_atom = atoms_unit_cell[a];
    
    new_atom[0] += L_unit_cell[0] * ix;
    new_atom[1] += L_unit_cell[1] * iy;
    new_atom[2] += L_unit_cell[2] * iz;
    
    atoms_.push_back(new_atom);
  }
}


vector<vector<double>> Crystal_Lattice::get_atoms_scaled(double dnn) const
{
  vector<vector<double>> atoms = atoms_;
  
  for (int i=0; i<atoms.size(); i++)
  {
    atoms[i][0] *= dnn;
    atoms[i][1] *= dnn;
    atoms[i][2] *= dnn;
  }
  
  return atoms;
}


vector<vector<double>> Crystal_Lattice::get_atoms_scaled(double Lx, double Ly, double Lz) const
{
  vector<vector<double>> atoms = atoms_;
  
  for (int i=0; i<atoms.size(); i++)
  {
    atoms[i][0] *= Lx/L_[0];
    atoms[i][1] *= Ly/L_[1];
    atoms[i][2] *= Lz/L_[2];
  }
  
  return atoms;
}


void Crystal_Lattice::export_to_xyz_file(string filename) const
{
  const int prec  = 8;
  const int ncols = 8+prec;
  
  ofstream of1(filename.c_str());
  
  of1 << atoms_.size() << endl;
  of1 << fixed << setprecision(prec);
  of1 << "Lattice=\"" << L_[0] << " 0.0" << " 0.0" << " 0.0 " << L_[1] << " 0.0" << " 0.0" << " 0.0 " << L_[2] << "\"" << endl;
  
  for (int i=0; i<atoms_.size(); i++)
  {
    of1 << "Atom"
        << setw(ncols) << atoms_[i][0]
        << setw(ncols) << atoms_[i][1]
        << setw(ncols) << atoms_[i][2]
        << endl;
  }
  
  of1.close();
}












