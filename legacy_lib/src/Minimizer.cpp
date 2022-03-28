#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <chrono>

using namespace std;

#include "Minimizer.h"
#include "DFT.h"
#include "myColor.h"

void Minimizer::run(long maxSteps)
{
  initialize();

  resume(maxSteps);
}

void Minimizer::resume(long maxSteps)
{
  do {
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    
    F_ = step();

    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
    elapsed_seconds_ = now-start;
    
    step_counter_++;

    double Ntotal = dft_->getNumberAtoms(0);
    double Volume = dft_->lattice().getVolume();

    f_abs_max_ = get_convergence_monitor();

    draw_after();

    if(f_abs_max_ < forceLimit_)
      {
	stringstream s;	
	s << "dF sufficiently small ... normal exit";
	reportMessage(s.str());
	break;
      }
    if(maxSteps > 0 && step_counter_ == maxSteps)
      {
	stringstream s;
	s << "maxSteps = " << maxSteps << " reached ... normal exit";
	reportMessage(s.str());
	break;
      }
    if(Ntotal/Volume < minDensity_)
      {
	stringstream s;
	s << "Density is " << Ntotal/Volume << " < minDensity_ = " << minDensity_ << " ... exiting" << endl;
	reportMessage(s.str());
	break;
      }
  } while(1);

}


void Minimizer::initialize()
{
  step_counter_ = 0;
  calls_ = 0;

  for(int Jspecies=0;Jspecies<dft_->getNumberOfSpecies();Jspecies++)
    dft_->getSpecies(Jspecies)->get_density_alias(x_[Jspecies]);    
}


DFT_Vec Minimizer::getDF(int Jspecies)
{
	DFT_Vec dF; dF.set(dft_->getDF(Jspecies));
	
	const Density& theDensity = dft_->getDensity(Jspecies);
	FreezingParameters fparams = theDensity.get_freezing_parameters();
	
	long pos = 0;
	for (int ix=0; ix<theDensity.Nx(); ix++)
	for (int iy=0; iy<theDensity.Ny(); iy++)
	for (int iz=0; iz<theDensity.Nz(); iz++)
	{
		bool is_in_selection = true;
		if (ix<fparams.Nx_min || ix>fparams.Nx_max) is_in_selection = false;
		if (iy<fparams.Ny_min || iy>fparams.Ny_max) is_in_selection = false;
		if (iz<fparams.Nz_min || iz>fparams.Nz_max) is_in_selection = false;
		
		bool freeze = false;
		if ( is_in_selection &&  fparams.freeze_inside_of_selection) freeze = true;
		if (!is_in_selection && !fparams.freeze_inside_of_selection) freeze = true;
		
		if (freeze) dF.set(pos, 0.0);
		pos++;
	}
	
	return dF;
}


double Minimizer::getDF_DX()
{
  calls_++;

  for(int Jspecies=0;Jspecies<x_.size();Jspecies++)
    dft_->getSpecies(Jspecies)->set_density_from_alias(x_[Jspecies]);        
  
  double F = 0;
  try {
    F = dft_->calculateFreeEnergyAndDerivatives(false);   
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  for(int Jspecies = 0; Jspecies<dft_->getNumberOfSpecies(); Jspecies++)
    dft_->getSpecies(Jspecies)->convert_to_alias_deriv(x_[Jspecies],dft_->getDF(Jspecies));    
  
  return F;
}

fireMinimizer2::fireMinimizer2(DFT *dft) :  Minimizer(dft)
{
  v_.resize(dft_->getNumberOfSpecies());
  for(auto &v: v_) v.resize(dft_->lattice().Ntot());

  dt_ = 1e-3; // my guess
  dt_max_ = 10*dt_;
    
  N_delay_ = 5;
  
  alpha_start_ = 0.1;
  f_dec_ = 0.5;
  f_inc_ = 1.1;
  f_alf_ = 0.99;  
}

void fireMinimizer2::initialize()
{
  //  fireMinimizer_Mu::initialize();

  Minimizer::initialize();
  
  it_    = 0;
  cut_   = 0;
  alpha_ = alpha_start_;

  N_P_positive_ = 0;
  N_P_negative_ = 0;

  vv_ = 1.0;

  dt_best_ = dt_;
  
  F_ = getDF_DX();
  
  for(auto &v: v_)
    v.zeros();

}


double fireMinimizer2::step()
{

  it_++;
  static double fold = F_;

  int begin_relax = 0;
  int end_relax   = dft_->getNumberOfSpecies();

  if(onlyRelax_ >= 0) {begin_relax = onlyRelax_; end_relax = begin_relax+1;}

  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  Summation P;  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    P += -v_[Jspecies].dotWith(getDF(Jspecies));

  //  if(F_ - fold > 1e-10) P = -1;
  fold = F_;

  int numSpecies = end_relax-begin_relax;
  vector<DFT_Vec> x_rem(numSpecies);
  vector<DFT_Vec> v_rem(numSpecies);
  vector<DFT_Vec> dF_rem(numSpecies);

  for(int Jspecies = begin_relax; Jspecies < end_relax; Jspecies++)
    {
      x_rem[Jspecies].set(x_[Jspecies]);
      v_rem[Jspecies].set(v_[Jspecies]);
      dF_rem[Jspecies].set(getDF(Jspecies));
    }

  
  if(P > 0)
    {
      N_P_positive_++;
      N_P_negative_ = 0;
      if(N_P_positive_ > N_delay_)
	{
	  dt_ = min(dt_*f_inc_,dt_max_);
	  alpha_ =alpha_*f_alf_;
	}
    } else {
    N_P_positive_ = 0;
    N_P_negative_++;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // This is the new, alternative stopping criterion: if there are too many steps in the wrong direction, hang up.
    // This needs to be communicated up the ladder so that it can be dealt with accordingly ... 
    if(N_P_negative_ > N_P_negative_max_) throw std::runtime_error("Cannot stop going uphill in Fire2");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Here, there is the possiblity of an initial delay of Ndelay_ steps before beginning the decrease of dt and of alpha. 
    if(!(initial_delay_ && it_ < N_delay_))
      {
	if(dt_*f_dec_ >= dt_min_)
	  dt_ *= f_dec_;
	alpha_ = alpha_start_;
      }
    reportMessage("Uphill motion stopped");
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
      {
	x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies], -0.5*dt_);
	v_[Jspecies].MultBy(0.1);
      }
    backtracks_++;
  }

  // integration step
  // Changed handling of backtracking 21/10/2019
  double rem = fmax_;
  try {
    SemiImplicitEuler(begin_relax, end_relax);
    vv_ = vv_*0.9 + 0.1*(fabs(F_ - fold)/dt_); // a measure of the rate at which the energy is changing
    dt_best_ = max(dt_,dt_best_);
  } catch(Eta_Too_Large_Exception &e) {
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
      {
	x_[Jspecies].set(x_rem[Jspecies]);
	v_[Jspecies].set(v_rem[Jspecies]);
	v_[Jspecies].MultBy(0.5);		
	dft_->setDF(Jspecies,dF_rem[Jspecies]);
      }
    dt_ /= 2;
    //    dt_max_ *= f_back_; // old method of control of dt_max_
    fmax_ = 1000; //rem;
    backtracks_++;
  }
  catch(...) {
    reportMessage("Unknown exception ...");
  }

  if(backtracks_ >= 10) // new control of dt_max_
    {
      dt_max_ = min(dt_best_, 0.9*dt_max_);
      dt_best_ = 0;
      backtracks_ = 0;
    };

  cout << "dt_best_ = " << dt_best_ << " backtracks_ = " << backtracks_ << endl;    
  // write a snapshot
  static int ic = 0;
  if(ic % 10 == 0)
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
      {
	stringstream s;
	s << "snapshot_s" << Jspecies << ".dat";
	string of = s.str();
	dft_->writeDensity(Jspecies,of);
      }
  ic++;
  
  return F_;
}

void fireMinimizer2::SemiImplicitEuler(int begin_relax, int end_relax)
{
  // update velocities and prepare for mixing
  // N.B.: df is a gradient, not a force
  double vnorm = 0.0;
  double fnorm = 0.0;
  long   cnorm = 0;
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      DFT_Vec df = getDF(Jspecies);      
      v_[Jspecies].IncrementBy_Scaled_Vector(df, -dt_);
	  
      double v = v_[Jspecies].euclidean_norm();
      double f = df.euclidean_norm();

      vnorm += v*v;
      fnorm += f*f;
      cnorm += df.size();      
    }
  rms_force_ = sqrt(fnorm/cnorm);
  vnorm_ = sqrt(vnorm/cnorm);
  
  /// Do mixing & update x
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      v_[Jspecies].MultBy(1-alpha_);      
      v_[Jspecies].IncrementBy_Scaled_Vector(getDF(Jspecies), -alpha_*sqrt(vnorm/fnorm));
      x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies],dt_);
    }  
  // recalculate forces with back-tracking, if necessary
  bool bSuccess = false;
  try{  
    F_ = getDF_DX(); // get new forces
    bSuccess = true;
  } catch (Eta_Too_Large_Exception &e) {
    reportMessage("Backtrack .. ");
    throw(e);
  }

  fnorm = 0.0;
  cnorm = 0;
  double new_fmax = 0;
  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      DFT_Vec df = getDF(Jspecies);
      double f = df.euclidean_norm();
      new_fmax = max(new_fmax, df.inf_norm()/dft_->getDensity(Jspecies).dV());
      fnorm += f*f;
      cnorm += df.size();
      //      cout << "Jspecies = " << Jspecies << " fmax_ = " << fmax_ << " df.inf_norm = " << df.inf_norm() << " fnorm = " << fnorm << " cnorm = " << cnorm << endl;
    }
  rms_force_ = sqrt(fnorm/cnorm);
  fmax_ = new_fmax;
}
