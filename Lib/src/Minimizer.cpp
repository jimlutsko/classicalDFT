#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

#include "Minimizer.h"
#include "DFT.h"

void Minimizer::run(long maxSteps)
{
  initialize();
  
  log_ << "Initialized ... removing old images ... " << endl;

  log_ << "Initial free energy = " << F_ << endl;
  log_ << "step_counter\tF\tF/N\tF/V\tf_abs_max\tN\tDensity\tCalls" << endl;


  //  int ret = system("rm *_image_*.png");
  int image_counter = 0;
  do {
    draw_before();

    F_ = step();

    step_counter_++;

    double Ntotal = dft_.getNumberAtoms(0);
    double Volume = dft_.lattice().getVolume();

    f_abs_max_ = get_convergence_monitor();

    log_.precision(12);
    log_ << step_counter_ 
	<< "\t" << F_ 
	<< "\t" << F_/Ntotal
	<< "\t" << F_/Volume
	<< "\t" << f_abs_max_
	<< "\t" << Ntotal
	<< "\t" << Ntotal/Volume
      	<< "\t" << calls_
	<< endl;
    if(std::isnan(f_abs_max_))
	log_  << "INF detected" << endl; //: min = " << dF_.min() << " max = " << dF_.max() << endl;

    draw_after();

    if(f_abs_max_ < forceLimit_)
      {
	log_ << "dF sufficiently small ... normal exit" << endl;
	break;
      }
    if(maxSteps > 0 && step_counter_ == maxSteps)
      {
	log_ << "maxSteps reached ... normal exit" << endl;
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

  double dV = dft_.lattice().dV();

  double ff = 0;
  double im = 0;

  int ix = 0;
  int iy = 0;
  int iz = 0;

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      auto &f = dft_.getDF(Jspecies);
      for(long i=0;i<f.size();i++)
	{
	  bool onBoundary = false;

	  if(bFrozenBoundary_)
	    {
	      iz++;
	      if(iz == dft_.lattice().Nz())
		{
		  iz = 0;
		  iy++;
		  if(iy == dft_.lattice().Ny())
		    {
		      iy = 0;
		      ix++;		  
		    }
		}
      
	      if(ix == 0 || iy == 0 || iz == 0)
		onBoundary = true;
	    }

      
	  double force = f.get(i);

	  if(!onBoundary)
	    f.set(i,2*force*x_[Jspecies].get(i));
	  else f.set(i,0.0);	
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

void fireMinimizer_Mu::draw_after()
{
  //  log_ << "After FIRE step " << step_counter_ << " F = " << F_ << " N = " << dft_.getNumberAtoms(0) << " calls = " << calls_ << " dt_ = " << dt_ << " alpha = " << alpha_ << endl;
}

int fireMinimizer_Mu::draw_during()
{
  log_ << "\t1D minimization F-mu*N = " << F_  << " N = " << dft_.getNumberAtoms(0) << endl;
  return 1;
}



double fireMinimizer2::step()
{
  it_++;

  int begin_relax = 0;
  int end_relax   = dft_.getNumberOfSpecies();

  if(onlyRelax_ >= 0) {begin_relax = onlyRelax_; end_relax = begin_relax+1;}

  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  double P = 0;  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    P += -v_[Jspecies].dotWith(dft_.getDF(Jspecies));

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
      {
	x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies],-0.5*dt_);
	v_[Jspecies].zeros(v_[Jspecies].size());
      }
  }

  // integration step
  
  bool blewUp = false;
  
  do {  
    try {
      blewUp = false;
      SemiImplicitEuler();
    } catch(Eta_Too_Large_Exception &e) {
      dt_ /= 2;
      dt_max_ /= 2;
      blewUp = true;
    }
  } while(blewUp);


  // write a snapshot
  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      stringstream s;
      s << "snapshot_s" << Jspecies << ".dat";
      string of = s.str();
      dft_.writeDensity(Jspecies,of);
    }
  
  return F_;
}

void fireMinimizer2::SemiImplicitEuler()
{
  int begin_relax = 0;
  int end_relax   = dft_.getNumberOfSpecies();

  if(onlyRelax_ >= 0) {begin_relax = onlyRelax_; end_relax = begin_relax+1;}

  // update velocities and prepare for mixing
  // N.B.: df is a gradient, not a force
  double vnorm = 0.0;
  double fnorm = 0.0;
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      DFT_Vec &df = dft_.getDF(Jspecies);
      v_[Jspecies].IncrementBy_Scaled_Vector(df, -dt_);    

      double v = v_[Jspecies].euclidean_norm();
      double f = df.euclidean_norm();

      vnorm += v*v;
      fnorm += f*f;
    }

  /// Do mixing
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      v_[Jspecies].MultBy(1-alpha_);      
      v_[Jspecies].IncrementBy_Scaled_Vector(dft_.getDF(Jspecies), -alpha_*sqrt(vnorm/fnorm));
    }

  //Update x
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)  
    x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies],dt_);           

  // recalculate forces with back-tracking, if necessary
  bool bSuccess = false;
  while(!bSuccess)
    {
      try{  
	F_ = getDF_DX();                          // then get new forces
	bSuccess = true;
      } catch (Eta_Too_Large_Exception &e) {
	for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)  
	  x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies],-0.5*dt_);
	dt_ /= 2;

	log_ << "Backtrack .. " << endl;
      }
    }
}
