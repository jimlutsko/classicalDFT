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
       << setw(20) << "dF/dt "
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
	 << setw(20) << vv_
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

  for(int Jspecies=0;Jspecies<dft_.getNumberOfSpecies();Jspecies++)
    {
      const Density& density = dft_.getDensity(Jspecies);
      x_[Jspecies].setAliasFromValues(density.getDensity());
    }
}


double Minimizer::getDF_DX()
{
  calls_++;

  for(int Jspecies=0;Jspecies<x_.size();Jspecies++)
    dft_.set_density_from_amplitude(Jspecies,x_[Jspecies]);

  double F = 0;
  try {
    F = dft_.calculateFreeEnergyAndDerivatives(false);   
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    dft_.getDF(Jspecies).alias_Jacobian(x_[Jspecies]);

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


fireMinimizer2::fireMinimizer2(DFT &dft, ostream &log) :  Minimizer(dft, log)
{
  v_.resize(dft_.getNumberOfSpecies());
  for(auto &v: v_) v.resize(dft_.lattice().Ntot());

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

  steps_since_backtrack_ = 0;
  
  F_ = getDF_DX();

  for(auto &v: v_)
    v.zeros();

}

double fireMinimizer2::step()
{

  //  if(steps_since_backtrack_ >= 1000) {log_ << "Increasing dt_max_" << endl; dt_max_ *= 10; steps_since_backtrack_ = 0;}
  
  it_++;
  static double fold = F_;

  int begin_relax = 0;
  int end_relax   = dft_.getNumberOfSpecies();

  if(onlyRelax_ >= 0) {begin_relax = onlyRelax_; end_relax = begin_relax+1;}

  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  Summation P;  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    P += -v_[Jspecies].dotWith(dft_.getDF(Jspecies));

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
    log_ << "Uphill motion stopped" << endl;
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
      {
	x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies], -0.5*dt_);
	//	v_[Jspecies].zeros();
	v_[Jspecies].MultBy(0.1);
	//v_[Jspecies].MultBy(0.5);
      }
  }

  // integration step
  // Changed handling of backtracking 21/10/2019  
  try {
    SemiImplicitEuler(begin_relax, end_relax);
    vv_ = vv_*0.9 + 0.1*(fabs(F_ - fold)/dt_); // a measure of the rate at which the energy is changing
    steps_since_backtrack_++;
  } catch(Eta_Too_Large_Exception &e) {
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
      {
	x_[Jspecies].set(x_rem[Jspecies]);
	v_[Jspecies].set(v_rem[Jspecies]);
	v_[Jspecies].MultBy(0.5);		
	dft_.setDF(Jspecies,dF_rem[Jspecies]);
      }
    //    dt_ /= 10; // This was previously 2
    dt_ /= 2;
    dt_max_ *= f_back_;
    steps_since_backtrack_ = 0;
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
  long   cnorm = 0;
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      DFT_Vec &df = dft_.getDF(Jspecies);      
      //      v_[Jspecies].IncrementBy_Scaled_Vector(df, -dt_);
      for(long i=0;i<v_[Jspecies].size();i++)
	{
	  double fac = -dt_;
	  //	  if(x_[Jspecies].get(i) < 0) fac *= (1-x_[Jspecies].get(i))*(1-x_[Jspecies].get(i)); 
	  //	  if(x_[Jspecies].get(i) >1) fac /= (x_[Jspecies].get(i) * x_[Jspecies].get(i));
	  v_[Jspecies].set(i, v_[Jspecies].get(i)+df.get(i)*fac);
	}
	  
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
      v_[Jspecies].IncrementBy_Scaled_Vector(dft_.getDF(Jspecies), -alpha_*sqrt(vnorm/fnorm));
      //      x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies],dt_);
      DFT_Vec &df = dft_.getDF(Jspecies);
      for(long i=0;i<x_[Jspecies].size();i++)
	{
	  double fac = dt_;
	  //	  if(x_[Jspecies].get(i) < 0) fac *= (1-x_[Jspecies].get(i)) * (1-x_[Jspecies].get(i)); 
	  //	  if(x_[Jspecies].get(i) >1) fac /= (x_[Jspecies].get(i) * x_[Jspecies].get(i)); 	  
	  x_[Jspecies].set(i, x_[Jspecies].get(i)+v_[Jspecies].get(i)*fac);
	}      
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

  fnorm = 0.0;
  cnorm = 0;
  fmax_ = 0;
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      DFT_Vec &df = dft_.getDF(Jspecies);
      double f = df.euclidean_norm();
      fmax_ = max(fmax_, df.inf_norm()/dft_.getDensity(Jspecies).dV());
      fnorm += f*f;
      cnorm += df.size();            
      if(threshold_ > 0)	
	for(long i=0;i<df.size();i++)
	  {

	    double fac = 1;
	    if(x_[Jspecies].get(i) < 0) fac *= 10;

	    //	    if(x_[Jspecies].get(i) < 0) fac *= (1-x_[Jspecies].get(i)); // * (1-x_[Jspecies].get(i)); // * (1-x_[Jspecies].get(i));
	    //	    if(x_[Jspecies].get(i) >1) fac /= (10*x_[Jspecies].get(i)); // * x_[Jspecies].get(i));
	    if(x_[Jspecies].get(i) > 0) fac /= exp(x_[Jspecies].get(i)); // * x_[Jspecies].get(i));
	    df.set(i, df.get(i) * fac);
	  }
    }
    rms_force_ = sqrt(fnorm/cnorm);
}
