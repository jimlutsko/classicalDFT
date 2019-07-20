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


void Minimizer::run(string& logfile, long maxSteps)
{
  initialize();
  
  cout << "Initialized ... removing old images ... " << endl;

  ofstream log(logfile.c_str(),ios::app);
  log << "#Initial free energy = " << F_ << endl;
  log << "#step_counter\tF\tF/N\tF/V\tf_abs_max\tN\tDensity\tCalls" << endl;
  log.close();

  int ret = system("rm image_*.png");
  int image_counter = 0;
  do {
    draw_before();

    F_ = step();

    step_counter_++;

    double Ntotal = dft_.getNumberAtoms(0);
    double Volume = dft_.lattice().getVolume();

    f_abs_max_ = get_convergence_monitor();

    ofstream log(logfile.c_str(), ios::app);
    log.precision(12);
    log << step_counter_ 
	<< "\t" << F_ 
	<< "\t" << F_/Ntotal
	<< "\t" << F_/Volume
	<< "\t" << f_abs_max_
	<< "\t" << Ntotal
	<< "\t" << Ntotal/Volume
      	<< "\t" << calls_
	<< endl;
    if(std::isnan(f_abs_max_))
      {
	log  << "INF detected" << endl; //: min = " << dF_.min() << " max = " << dF_.max() << endl;
	cout << "INF detected" << endl; //: min = " << dF_.min() << " max = " << dF_.max() << endl;
      }
    log.close();

    draw_after();

    int oldprecision = cout.precision(12);
    cout << "F = " << F_ << " f_abs_max_ = " << f_abs_max_ << endl;
    cout.precision(oldprecision);

    if(step_counter_%20 == 1)
      {
	stringstream ss;
	std::ios  state(NULL);
	state.copyfmt(ss);
	
	ss << "image_";
	ss << setfill ('0') << std::setw(8);
	ss << image_counter;
	ss.copyfmt(state);
	ss <<".png";
	finish(ss.str().c_str());
	image_counter++;
      } else {
      	stringstream ss;
	ss << "image_current" <<".png";
	finish(ss.str().c_str());
    }
    /*
    if(step_counter_%1000 == 1)
      {
	stringstream ss;
	ss << "snapshot_" << step_counter_ <<".dat";
	string s = ss.str();
	density_.writeDensity(s); //s.str());
      }
    */

    if(f_abs_max_ < forceLimit_)
      {
	cout << "dF sufficiently small ... normal exit" << endl;
	break;
      }
    if(maxSteps > 0 && step_counter_ == maxSteps)
      {
	cout << "maxSteps reached ... normal exit" << endl;
	break;
      }
  } while(1);

  finish("final.png");
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
  //  dF_.multBy(2);

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


void Minimizer_Fixed_N::initialize()
{
  Minimizer::initialize();
}

double Minimizer_Fixed_N::getDF_DX()
{
  calls_++;

  double dV = dft_.lattice().dV();
  vector<double> sum;
  vector<double> a;

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      double N0 = x_[Jspecies].size()*SMALL_VALUE*dV;
      sum.push_back(x_[Jspecies].dotWith(x_[Jspecies])); 

      a.push_back((N_fixed_target_-N0)/(sum.back()*dV));

      DFT_Vec y(x_[Jspecies]);
      y.multBy(sqrt(a.back()));

      dft_.set_density_from_amplitude(Jspecies, y);
    }

  double F = 0;

  try {
    F = dft_.calculateFreeEnergyAndDerivatives(false);   
    cout.precision(12);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  double ff = 0;

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      DFT_Vec &df = dft_.getDF(Jspecies);
      
      mu_eff_ = 0.0;
      for(long i=0;i<df.size();i++)
	mu_eff_ += df.get(i)*x_[Jspecies].get(i)*x_[Jspecies].get(i);
      mu_eff_ *= 1.0/sum[Jspecies];

      for(int i=0;i<df.size();i++)
	{
	  double force = df.get(i)-mu_eff_;
	  df.set(i,2*a[Jspecies]*force*x_[Jspecies].get(i));  
	}

    }
  return F;
}


//static double mu1;

void fireMinimizer::initialize()
{
  Minimizer_Fixed_N::initialize();
  
  it_ = 0;
  cut_ = 0;

  alpha_ = alpha_start_;
  
  F_ = getDF_DX();

  for(auto &v: v_)
    v.zeros(v.size());
}


void fireMinimizer::verlet()
{
  vector<DFT_Vec> y;
  for(auto &x: x_) y.push_back(x);
    
  try{
    for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
      {
	DFT_Vec &df = dft_.getDF(Jspecies);
	
	x_[Jspecies].Increment_And_Scale(v_[Jspecies],dt_);// this gives x_(t+dt)
	x_[Jspecies].Increment_And_Scale(df,-0.5*dt_*dt_); // dF=dV/dx does not have minus sign ... 
	v_[Jspecies].Increment_And_Scale(df,-0.5*dt_);     // now do half the v update
      }
    F_ = getDF_DX();                          // then get new forces    
    for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
      v_[Jspecies].Increment_And_Scale(dft_.getDF(Jspecies), -0.5*dt_);    // and finish velocity update
  } catch (Eta_Too_Large_Exception &e) {
    for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
      {
	x_[Jspecies].set(y[Jspecies]);
	v_[Jspecies].zeros();
      }
    cout << "Backtrack .. " << endl;
    throw e;
  }
}


double fireMinimizer::step()
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
      stringstream s;
      s << "snapshot_s" << Jspecies << ".dat";
      string of = s.str();
      dft_.writeDensity(Jspecies,of);
    }  
  
  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  double p = 0;
  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      p += -v_[Jspecies].dotWith(dft_.getDF(Jspecies));
      v_[Jspecies].multBy(1-alpha_);
    }

  double vnorm = 0;
  double fnorm = 0;

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      double v = v_[Jspecies].euclidean_norm();
      double f = dft_.getDF(Jspecies).euclidean_norm();
      vnorm += v*v;
      fnorm += f*f;
    }

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)  
    v_[Jspecies].Increment_And_Scale(dft_.getDF(Jspecies), -alpha_*sqrt(vnorm/fnorm));

  if(p < 0)
    {
      cout << "\tp < 0 : reset" << endl;
      for(auto &v: v_)
	v.zeros(v.size());
      cut_ = it_;
      dt_ *= f_dec_;
      alpha_ = alpha_start_;
    } else if(it_-cut_>N_min_) {
    dt_ = min(dt_*f_inc_,dt_max_);
    alpha_ =alpha_*f_alf_;
  }
  return F_;
}

void fireMinimizer::draw_after()
{
  cout << "After FIRE step " << step_counter_ << " F = " << F_ << " N = " << dft_.getNumberAtoms(0) << " calls = " << calls_ << " dt_ = " << dt_ << " alpha = " << alpha_ << endl;
}

int fireMinimizer::draw_during()
{
  cout << "\t1D minimization F = " << F_  << " mu_eff = " << mu_eff_ << " F-mu*N = " << F_-mu_eff_*N_fixed_target_ << " N = " << dft_.getNumberAtoms(0) << " Ntarget = " << N_fixed_target_ << endl;
  return 1;
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
	DFT_Vec &df = dft_.getDF(Jspecies);
	
	x_[Jspecies].Increment_And_Scale(v_[Jspecies],dt_);           // this gives x_(t+dt)
	x_[Jspecies].Increment_And_Scale(df,-0.5*dt_*dt_); // dF=dV/dx does not have minus sign ... 
	v_[Jspecies].Increment_And_Scale(df, -0.5*dt_);    // now do half the v update
      }
    F_ = getDF_DX();                          // then get new forces    
    for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
      v_[Jspecies].Increment_And_Scale(dft_.getDF(Jspecies), -0.5*dt_);    // and finish velocity update
  } catch (Eta_Too_Large_Exception &e) {
    for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
      {
	x_[Jspecies].set(y[Jspecies]);
	v_[Jspecies].zeros(v_[Jspecies].size());
      }
    cout << "Backtrack .. " << endl;
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
      stringstream s;
      s << "snapshot_s" << Jspecies << ".dat";
      string of = s.str();
      dft_.writeDensity(Jspecies,of);
    }
  
  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  double p = 0;
  
  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      p += -v_[Jspecies].dotWith(dft_.getDF(Jspecies));
      v_[Jspecies].multBy(1-alpha_);
    }

  double vnorm = 0;
  double fnorm = 0;

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)
    {
      double v = v_[Jspecies].euclidean_norm();
      double f = dft_.getDF(Jspecies).euclidean_norm();
      vnorm += v*v;
      fnorm += f*f;
    }

  for(int Jspecies = 0; Jspecies<dft_.getNumberOfSpecies(); Jspecies++)  
    v_[Jspecies].Increment_And_Scale(dft_.getDF(Jspecies), -alpha_*sqrt(vnorm/fnorm));

  if(p < 0)
    {
      cout << "\tp < 0 : reset" << endl;
      for(auto &v: v_)
	v.zeros(v.size());
      cut_ = it_;
      dt_ *= f_dec_;
      alpha_ = alpha_start_;
    } else if(it_-cut_>N_min_) {
    dt_ = min(dt_*f_inc_,dt_max_);
    alpha_ =alpha_*f_alf_;
  }
  return F_;
}

void fireMinimizer_Mu::draw_after()
{
  cout << "After FIRE step " << step_counter_ << " F = " << F_ << " N = " << dft_.getNumberAtoms(0) << " calls = " << calls_ << " dt_ = " << dt_ << " alpha = " << alpha_ << endl;
}

int fireMinimizer_Mu::draw_during()
{
  cout << "\t1D minimization F-mu*N = " << F_  << " N = " << dft_.getNumberAtoms(0) << endl;
  return 1;
}
