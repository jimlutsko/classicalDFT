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
  
  log_ << "Initialized ... removing old images ... " << endl;

  log_ << "Initial free energy = " << F_ << endl;
  log_ << setw(12) <<"step_counter"
       << setw(20) << myColor::YELLOW << "F    "  << myColor::RESET
       << setw(20) << "F/N    "
       << setw(20) << "F/V    "
       << setw(20) << myColor::RED << "f_abs_max    "  << myColor::RESET
       << setw(20) << "N    "
       << setw(20) << "Density    "
       << setw(10) << myColor::RED << "Time(sec)"  << myColor::RESET
       << endl;
  
  int image_counter = 0;

  do {
    draw_before();

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    
    F_ = step();

    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = now-start;
    
    step_counter_++;

    double Ntotal = dft_.getNumberAtoms(0);
    double Volume = dft_.lattice().getVolume();

    f_abs_max_ = get_convergence_monitor();

    log_.precision(12);    
    log_ << setw(12) << step_counter_
	 << setw(20) << myColor::YELLOW << F_  << myColor::RESET
	 << setw(20) << F_/Ntotal
	 << setw(20) << F_/Volume
	 << setw(20) << myColor::RED << f_abs_max_  << myColor::RESET
	 << setw(20) << Ntotal
	 << setw(20) << Ntotal/Volume;
    log_.precision(4);        
    log_ << setw(10) << myColor::RED << elapsed_seconds.count() << myColor::RESET
	 << endl;

    
    if(std::isnan(f_abs_max_))
	log_  << "INF detected" << endl; 

    draw_after();

    if(f_abs_max_ < forceLimit_)
      {
	log_ << "dF sufficiently small ... normal exit" << endl;
	break;
      }
    if(maxSteps > 0 && step_counter_ == maxSteps)
      {
	log_ << "maxSteps = " << maxSteps << " reached ... normal exit" << endl;
	break;
      }
  } while(1);

}


void Minimizer::initialize()
{
  step_counter_ = 0;
  calls_ = 0;

  for(int species=0;species<dft_.getNumberOfSpecies();species++)
    {
      const Density& density = dft_.getDensity(species);
      for(long i=0;i<density.Ntot();i++)
	x_[species].set(i,sqrt(max(0.0, density.getDensity(i) - SMALL_VALUE)));      
    }
}


double Minimizer::getDF_DX()
{
  calls_++;

  for(int i=0;i<x_.size();i++)
    dft_.set_density_from_amplitude(i,x_[i]);

  double F = 0;
  try {
    F = dft_.calculateFreeEnergyAndDerivatives(false);   
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  //  dF_.Schur(x_,dF_);
  //  dF_.MultBy(2);


  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      auto &f = dft_.getDF(Jspecies);
      for(long i=0;i<f.size();i++)
	{
	  double force = f.get(i);
	  f.set(i,2*force*x_[Jspecies].get(i));
	}
    }
  return F;
}

void fireMinimizer_Mu::initialize()
{
  Minimizer::initialize();
  
  it_ = 0;
  cut_ = 0;

  alpha_ = alpha_start_;
  
  F_ = getDF_DX();

  for(auto &v: v_)
    v.zeros(v.size());
}


void fireMinimizer_Mu::verlet()
{
  vector<DFT_Vec> y;
  for(auto &x: x_) y.push_back(x);
    
  try{
    for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
      {
	if(onlyRelax_ >= 0)
	  if(Jspecies != onlyRelax_)
	    continue;
	DFT_Vec &df = dft_.getDF(Jspecies);
	
	x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies],dt_);           // this gives x_(t+dt)
	x_[Jspecies].IncrementBy_Scaled_Vector(df,-0.5*dt_*dt_); // dF=dV/dx does not have minus sign ... 
	v_[Jspecies].IncrementBy_Scaled_Vector(df, -0.5*dt_);    // now do half the v update
      }
    F_ = getDF_DX();                          // then get new forces    
    for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
      	if(onlyRelax_ >= 0)
	  if(Jspecies != onlyRelax_)
	    continue;
	  else v_[Jspecies].IncrementBy_Scaled_Vector(dft_.getDF(Jspecies), -0.5*dt_);    // and finish velocity update
  } catch (Eta_Too_Large_Exception &e) {
    for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
      {
	if(onlyRelax_ >= 0)
	  if(Jspecies != onlyRelax_)
	    continue;
	
	x_[Jspecies].set(y[Jspecies]);
	v_[Jspecies].zeros(v_[Jspecies].size());
      }
   log_ << "Backtrack .. " << endl;
    throw e;
  }
}


double fireMinimizer_Mu::step()
{
  it_++;

  bool blewUp = false;
  
  do {  
    try {
      blewUp = false;
      verlet();
    } catch(Eta_Too_Large_Exception &e) {
      dt_ /= 2;
      dt_max_ /= 2;
      blewUp = true;
    }
  } while(blewUp);


  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      if(onlyRelax_ >= 0)
	if(Jspecies != onlyRelax_)
	  continue;
	
      stringstream s;
      s << "snapshot_s" << Jspecies << ".dat";
      string of = s.str();
      dft_.writeDensity(Jspecies,of);
    }
  
  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  double p = 0;
  
  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      if(onlyRelax_ >= 0)
	if(Jspecies != onlyRelax_)
	  continue;
      
      p += -v_[Jspecies].dotWith(dft_.getDF(Jspecies));
      v_[Jspecies].MultBy(1-alpha_);
    }

  double vnorm = 0;
  double fnorm = 0;

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      if(onlyRelax_ >= 0)
	if(Jspecies != onlyRelax_)
	  continue;
	
      double v = v_[Jspecies].euclidean_norm();
      double f = dft_.getDF(Jspecies).euclidean_norm();
      vnorm += v*v;
      fnorm += f*f;
    }

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)  
    v_[Jspecies].IncrementBy_Scaled_Vector(dft_.getDF(Jspecies), -alpha_*sqrt(vnorm/fnorm));

  if(p < 0)
    {
      for(auto &v: v_)
	v.zeros(v.size());
      cut_ = it_;
      dt_ *= f_dec_;
      alpha_ = alpha_start_;
    } else if(it_-cut_>N_delay_) {
    dt_ = min(dt_*f_inc_,dt_max_);
    alpha_ =alpha_*f_alf_;
  }
  return F_;
}

void fireMinimizer_Mu::draw_after() {}

int fireMinimizer_Mu::draw_during()
{
  log_ << "\t1D minimization F-mu*N = " << F_  << " N = " << dft_.getNumberAtoms(0) << endl;
  return 1;
}



double fireMinimizer2::step()
{
  it_++;
  static double fold = F_;

  int begin_relax = 0;
  int end_relax   = dft_.getNumberOfSpecies();

  if(onlyRelax_ >= 0) {begin_relax = onlyRelax_; end_relax = begin_relax+1;}

  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  Summation P;  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    P += -v_[Jspecies].dotWith(dft_.getDF(Jspecies));

  if(F_ - fold > 1e-10) P = -1;
  fold = F_;

  int numSpecies = end_relax-begin_relax;
  vector<DFT_Vec> x_rem(numSpecies);
  vector<DFT_Vec> v_rem(numSpecies);
  vector<DFT_Vec> dF_rem(numSpecies);

  for(int Jspecies = begin_relax; Jspecies < end_relax; Jspecies++)
    {
      x_rem[Jspecies].set(x_[Jspecies]);
      v_rem[Jspecies].set(v_[Jspecies]);
      dF_rem[Jspecies].set(dft_.getDF(Jspecies));
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

    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
	v_[Jspecies].zeros(v_[Jspecies].size());
  }

  // integration step
  // Changed handling of backtracking 21/10/2019  
  try {
    SemiImplicitEuler(begin_relax, end_relax);
  } catch(Eta_Too_Large_Exception &e) {
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
      {
	x_[Jspecies].set(x_rem[Jspecies]);
	v_[Jspecies].set(v_rem[Jspecies]);
	dft_.setDF(Jspecies,dF_rem[Jspecies]);
      }
    dt_ /= 10; // This was previously 2
  }

  // write a snapshot
  static int ic = 0;
  if(ic % 10 == 0)
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      stringstream s;
      s << "snapshot_s" << Jspecies << ".dat";
      string of = s.str();
      dft_.writeDensity(Jspecies,of);
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
  long   count = 0;
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      DFT_Vec &df = dft_.getDF(Jspecies);
      count += df.size();
      v_[Jspecies].IncrementBy_Scaled_Vector(df, -dt_);
	  
      double v = v_[Jspecies].euclidean_norm();
      double f = df.euclidean_norm();

      vnorm += v*v;
      fnorm += f*f;
    }
  rms_force_ = sqrt(fnorm/count);
  
  /// Do mixing
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      v_[Jspecies].MultBy(1-alpha_);      
      v_[Jspecies].IncrementBy_Scaled_Vector(dft_.getDF(Jspecies), -alpha_*sqrt(vnorm/fnorm));
    }
  
  //Update x
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      long Ntot = x_[Jspecies].size();
      int chunk = Ntot/20;
      long i;
#pragma omp parallel for		       \
            shared( chunk, Jspecies, v_)	\
            private(i)				\
            schedule(static,chunk)				
	for(i = 0; i<Ntot; i++)
	  x_[Jspecies].set(i, x_[Jspecies].get(i) + v_[Jspecies].get(i)*dt_/max(fudge_,fabs(x_[Jspecies].get(i)))); // This weighting reduces backtracking ...
    }
  
  // recalculate forces with back-tracking, if necessary
  bool bSuccess = false;
  try{  
    F_ = getDF_DX(); // get new forces
    bSuccess = true;
  } catch (Eta_Too_Large_Exception &e) {
    log_ << "Backtrack .. " << endl;
    throw(e);
  }
}
