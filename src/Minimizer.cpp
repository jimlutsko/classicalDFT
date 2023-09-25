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


extern int counter_new_fft_plans;
extern int counter_fft_transforms;


Minimizer::Minimizer(DFT *dft) : dft_(dft), forceLimit_(0.1)
{
  x_.resize(dft->getNumberOfSpecies());
  
  for(auto &x: x_) x.resize(dft_->get_lattice().Ntot());

  for(int Jspecies=0;Jspecies<dft_->getNumberOfSpecies();Jspecies++)
    dft_->getSpecies(Jspecies)->get_density_alias(x_[Jspecies]);          
}

void Minimizer::reset()
{
  for(int Jspecies=0;Jspecies<dft_->getNumberOfSpecies();Jspecies++)
    dft_->getSpecies(Jspecies)->get_density_alias(x_[Jspecies]);        
}

void Minimizer::run(long maxSteps)
{
  step_counter_ = 0;
  calls_ = 0;

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
    double Volume = dft_->get_lattice().getVolume();

    f_abs_max_ = get_convergence_monitor();

    draw_after();

    if(should_stop()) break;
    
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

  cleanup();
}

bool Minimizer::should_stop() const
{
  bool stop = false;
  
  if(f_abs_max_ < forceLimit_)
  {
    stringstream s;	
    s << "Finished: convergence monitor = " << f_abs_max_ << " < " << forceLimit_ << " = forceLimit_ and so is  sufficiently small ... normal exit";
    cout << s.str() << endl;
    reportMessage(s.str());
    stop = true;
  }
  
  return stop;
}

void Minimizer::reset_fixed_directions()
{
  fixed_directions_.clear();
}

void Minimizer::set_fixed_direction(const DFT_Vec& fixed, bool already_using_density_alias)
{
  fixed_directions_.clear();
  add_fixed_direction(fixed, already_using_density_alias);
  /*
  fixed_directions_.push_back(fixed);
  fixed_directions_[0].normalise();
  
  if (!already_using_density_alias)
  {
    // Convert to alias space as we later orthogonalise the forces/velocity in alias space
    dft_->getSpecies(0)->convert_to_alias_increment(x_[0],fixed_directions_[0]);
    fixed_directions_[0].normalise();
  }
  
  F_ = getDF_DX();*/
}

void Minimizer::add_fixed_direction(const DFT_Vec& fixed, bool already_using_density_alias)
{
  fixed_directions_.push_back(fixed);
  fixed_directions_.back().normalise();
  
  if (!already_using_density_alias)
  {
    // Convert to alias space as we later orthogonalise the forces/velocity in alias space
    dft_->getSpecies(0)->convert_to_alias_increment(x_[0],fixed_directions_.back());
    fixed_directions_.back().normalise();
  }
  
  F_ = getDF_DX();
}

double Minimizer::getDF_DX()
{
  dft_->set_densities_from_aliases(x_);
  
  double F = 0;
  try {F = get_energy_and_forces(); calls_++;} 
  catch( Eta_Too_Large_Exception &e) {throw e;}

  //  if(use_squared_forces_) dft_->dF_to_H_dot_dF();
  
  dft_->convert_dF_to_alias_derivs(x_);
  
  // Project the forces into the subspace orthogonal to "fixed_direction_",
  // or revert the sign of the normal force if "flip_forces_along_fixed_direction_"
  // is set to true.
  
  for (int i=0; i<fixed_directions_.size(); i++)
  if(fixed_directions_[i].size() == dft_->getDF(0).size())
  {
    for(int s=0;s<dft_->getNumberOfSpecies();s++)
    {
      DFT_Vec& df = dft_->getDF(s);
      
      double fac = 1; if (flip_forces_along_fixed_direction_) fac = 2;
      df.IncrementBy_Scaled_Vector(fixed_directions_[i], -fac*fixed_directions_[i].dotWith(df));
    }
  }
  
  return F;
}

fireMinimizer2::fireMinimizer2(DFT *dft) :  Minimizer(dft)
{
  v_.resize(dft_->getNumberOfSpecies());
  for(auto &v: v_) v.resize(dft_->get_lattice().Ntot());

  dt_ = 1e-3; // my guess
  dt_max_ = 10*dt_;
    
  N_delay_ = 5;
  
  alpha_start_ = 0.1;
  f_dec_ = 0.5;
  f_inc_ = 1.1;
  f_alf_ = 0.99;

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

void fireMinimizer2::reset()
{
  Minimizer::reset();
  dt_best_ = dt_;
  
  F_ = getDF_DX();
  
  for(auto &v: v_)
    v.zeros();  
}

double fireMinimizer2::step()
{
  it_++;
  static double fold = F_;
  
  cout << "Before new step in minimizer: counter_new_fft_plans = " << counter_new_fft_plans << " counter_fft_transforms = " << counter_fft_transforms << endl;

  int begin_relax = 0;
  int end_relax   = dft_->getNumberOfSpecies();

  if(onlyRelax_ >= 0) {begin_relax = onlyRelax_; end_relax = begin_relax+1;}

  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  Summation P;  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    P += -v_[Jspecies].dotWith(dft_->getDF(Jspecies));

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
    dF_rem[Jspecies].set(dft_->getDF(Jspecies));
  }
  
  //long Ntot = dft_->getDF(0).size();
  //cout << setprecision(8);
  //cout << "\t-dF*v/|v||dF| = " << setw(16) << P/vnorm_/rms_force_/Ntot
  //     << "  -dF*v/|dF| = "    << setw(16) << P/rms_force_/sqrt(Ntot)
  //     << "  -dF*v/|v| = "     << setw(16) << P/vnorm_/sqrt(Ntot)
  //     << endl;
  
  // 15/12/2022: Changed so as to ALWAYS accept first iteration since v_ should be zero ...
  cout << scientific << "P = " << P << fixed << "    s0: v*F/|v||F| = " << v_[0].dotWith(dft_->getDF(0))/v_[0].euclidean_norm()/dft_->getDF(0).euclidean_norm() << endl;
  if(P > 0 || it_ == 1)
    {
      N_P_positive_++;
      N_P_negative_ = 0;
      
      if(N_P_positive_ > N_delay_)
	{
	  dt_ = min(dt_*f_inc_,dt_max_);
	  alpha_ =alpha_*f_alf_;
	}
    } else  {
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
	if(dt_*f_dec_ >= dt_min_) dt_ *= f_dec_;
	alpha_ = alpha_start_;
      }
    
    reportMessage("Uphill motion stopped");
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
      {
	x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies], -0.5*dt_);
	//v_[Jspecies].MultBy(0.1); vnorm_ *= 0.1;
	//	v_[Jspecies].MultBy(0.01); vnorm_ *= 0.01;
	v_[Jspecies].MultBy(0.0); vnorm_ *= 0.0;	
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
  }  catch(Eta_Too_Large_Exception &e) {
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
      {
	x_[Jspecies].set(x_rem[Jspecies]);
	v_[Jspecies].set(v_rem[Jspecies]);
	v_[Jspecies].MultBy(0.5); vnorm_ *= 0.5;
	dft_->setDF(Jspecies,dF_rem[Jspecies]);
      }
    dt_ /= 2;
    //    dt_max_ *= f_back_; // old method of control of dt_max_
    fmax_ = 1000; //rem;
    backtracks_++;
  }  catch(...)  {
    reportMessage("Unknown exception ...");
  }
  
  if(backtracks_ >= 10) // new control of dt_max_
  {
    stringstream s;
    s << "Changinging dt_max_ from " << dt_max_ << " to " << min(dt_best_, 0.9*dt_max_) << endl;
    reportMessage(s.str());
    dt_max_ = min(dt_best_, 0.9*dt_max_);
    dt_best_ = 0;
    backtracks_ = 0;
  };

  //  cout << "dt_best_ = " << dt_best_ << " backtracks_ = " << backtracks_ << endl;    
  // write a snapshot
  static int ic = 0;
  if(ic % 10 == 0)
  {
    for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
    {
      stringstream s;
      s << "snapshot_s" << Jspecies << ".dat";
      string of = s.str();
      dft_->writeDensity(Jspecies,of);

      stringstream s1;
      s1 << "snapshot_s" << Jspecies << ".vtk";
      string of1 = s1.str();
      dft_->write_density_vtk(Jspecies,of1);      
    }
  }
  ic++;
  
  return F_;
}

void fireMinimizer2::SemiImplicitEuler(int begin_relax, int end_relax)
{
  // Update velocities and prepare for mixing
  // N.B.: df is a gradient, not a force
  double vnorm = 0.0;
  double fnorm = 0.0;
  long   cnorm = 0;
  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
  {
    DFT_Vec &df = dft_->getDF(Jspecies);      
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
    double v = v_[Jspecies].euclidean_norm();
    double f = dft_->getDF(Jspecies).euclidean_norm();
    
    v_[Jspecies].MultBy(1-alpha_);      
    v_[Jspecies].IncrementBy_Scaled_Vector(dft_->getDF(Jspecies), -alpha_*v/f);
    x_[Jspecies].IncrementBy_Scaled_Vector(v_[Jspecies],dt_);
  }
  
  // Recalculate forces with back-tracking, if necessary
  bool bSuccess = false;
  try {  
    F_ = getDF_DX(); // get new forces
    bSuccess = true;
  } catch (Eta_Too_Large_Exception &e)  {
    reportMessage("Backtrack .. ");
    throw(e);
  }

  vnorm = 0.0;
  fnorm = 0.0;
  cnorm = 0;
  double new_fmax = 0;
  
  for(int Jspecies = begin_relax; Jspecies<end_relax; Jspecies++)
  {
    DFT_Vec &df = dft_->getDF(Jspecies);
    
    new_fmax = max(new_fmax, df.inf_norm()/dft_->getDensity(Jspecies).dV());
    double f = df.euclidean_norm();
    double v = v_[Jspecies].euclidean_norm();
    
    fnorm += f*f;
    vnorm += v*v;
    cnorm += df.size();
  }
  
  rms_force_ = sqrt(fnorm/cnorm);
  vnorm_ = sqrt(vnorm/cnorm);
  fmax_ = new_fmax;
}
