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

    double Ntotal = density_.getNumberAtoms();
    double Volume = density_.getVolume();

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

  for(long i=0;i<density_.Ntot();i++)
    x_.set(i,sqrt(max(0.0, density_.getDensity(i) - SMALL_VALUE)));
}


double Minimizer::getDF_DX()
{
  calls_++;

  density_.set_density_from_amplitude(x_);

  double F = 0;

  try {
    F = dft_.calculateFreeEnergyAndDerivatives(density_,mu_, dF_,false);   
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  //  dF_.Schur(x_,dF_);
  //  dF_.multBy(2);

  double dV = density_.dV();

  double ff = 0;
  double im = 0;
  double dm = 0;

  int ix = 0;
  int iy = 0;
  int iz = 0;
  
  for(int i=0;i<dF_.size();i++)
    {
      bool onBoundary = false;
      
      if(bFrozenBoundary_)
	{
	  iz++;
	  if(iz == density_.Nz())
	    {
	      iz = 0;
	      iy++;
	      if(iy == density_.Ny())
		{
		  iy = 0;
		  ix++;		  
		}
	    }
      
	  if(ix == 0 || iy == 0 || iz == 0)
	    onBoundary = true;
	}

      
      double force = dF_.get(i);

      if(!onBoundary)
	{
	  double density = density_.getDensity(i);
	  double density_predicted = exp(log(density)-force/dV);
	  if(fabs(density-density_predicted) > ff) 
	    { ff = fabs(density-density_predicted); im = i;}
	  dF_.set(i,2*force*x_.get(i));
	} else dF_.set(i,0.0);	
    }
  err_ = ff;

  return F;
}


void Minimizer_Fixed_N::initialize()
{
  Minimizer::initialize();

  for(long i=0;i<density_.Ntot();i++)
    x_.set(i,sqrt(max(0.0, density_.getDensity(i) - SMALL_VALUE)));
}

double Minimizer_Fixed_N::getDF_DX()
{
  calls_++;

  double dV = density_.dV();
  double N0 = x_.size()*SMALL_VALUE*dV;
  double sum = x_.dotWith(x_); 

  double a = (N_fixed_target_-N0)/(sum*dV);

  DFT_Vec y(x_);
  y.multBy(sqrt(a));

  density_.set_density_from_amplitude(y);

  double F = 0;

  try {
    F = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF_,false);   
    cout.precision(12);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  mu_eff_ = 0.0;
  for(long i=0;i<dF_.size();i++)
    mu_eff_ += dF_.get(i)*x_.get(i)*x_.get(i);
  mu_eff_ *= 1.0/sum;

  double ff = 0;
  double im = 0;
  double dm = 0;
  for(int i=0;i<dF_.size();i++)
    {
      double force = dF_.get(i)-mu_eff_;

      double density = density_.getDensity(i);
      double density_predicted = exp(log(density)-force/dV);
      if(fabs(density-density_predicted) > ff) 
	{ ff = fabs(density-density_predicted); im = i;}
      dF_.set(i,2*a*force*x_.get(i));  
    }

  err_ = ff;

  return F;
}


void ConjugateGradients2::initialize()
{
  Minimizer::initialize();

  for(long i=0;i<density_.Ntot();i++)
    x_.set(i,sqrt(max(0.0, density_.getDensity(i) - SMALL_VALUE)));

  F_ = getDF_DX();


  cout << "Initial value of F = " << F_ << endl;

  r_.set(dF_); r_.multBy(-1); // = -dF_;
  d_.set(r_); // = r_;

  delta_new_ = r_.dotWith(d_);
  delta_0_ = delta_new_;
}

double ConjugateGradients2::getDF_DX()
{
  calls_++;

  density_.set_density_from_amplitude(x_);

  double F = 0;

  try {
    F = dft_.calculateFreeEnergyAndDerivatives(density_,mu_, dF_,false);   
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  //  dF_.Schur(x_,dF_);
  //  dF_.multBy(2);

  double dV = density_.dV();

  double ff = 0;
  double im = 0;
  double dm = 0;
  for(int i=0;i<dF_.size();i++)
    {
      double force = dF_.get(i);

      double density = density_.getDensity(i);
      double density_predicted = exp(log(density)-force/dV);
      if(fabs(density-density_predicted) > ff) 
	{ ff = fabs(density-density_predicted); im = i;}
      dF_.set(i,2*force*x_.get(i));
    }
  err_ = ff;


  cout << "Density = " << density_.getDensity(density_.Nx()/2, (1.5/4)*density_.Ny(), density_.Nz()/2) << endl;
  cout << "Force   = " << dF_.get(density_.pos(density_.Nx()/2, (1.5/4)*density_.Ny(), density_.Nz()/2)) << endl;

    cout << "Density = " << density_.getDensity(0, (1.5/4)*density_.Ny(), density_.Nz()/2) << endl;
  cout << "Force   = " << dF_.get(density_.pos(0, (1.5/4)*density_.Ny(), density_.Nz()/2)) << endl;

  return F;
}

void ConjugateGradients2::draw_after() 
{    
    cout << "After conj grad step " << step_counter_ 
	 << " F-mu*N = " << F_ 
	 << " and N = " << density_.getNumberAtoms() 
	 << endl;
}

int ConjugateGradients2::draw_during()
{
  double N = density_.getNumberAtoms();
  cout << "\t 1D minimization F-mu*N = " << F_ << " N = " << N << endl;
}

static double mu1;


/************************************************
The secant search is based on making a linear approximation:
f(x) = f(x0) + ((x-x0)/(x1-x0))*(f(x1)-f(x0))
giving the estimate for the zero as 

x = x0 - (x1-x0)*f(x0)/(f(x1)-f(x0))

Normally we would need
y0 = f(x)
x0 = x

x1 = x0 + alpha
y1 = f(x1)

x = x0 - (x1 - x0)*y0/(y1-y0)

Note that this can be written as
x-x1 = x0-x1 - (x1-x0)*y0/(y1-y0) = y1*(x1-x0)/(y0-y1)

so what we actually do is:
y0=f(x)

x += alpha
y1 = f(x)

alpha ==> alpha*(y1-y0)/(y0-y1)
***************************************************************/

double ConjugateGradients2::tryStep(double alf, double &dalf, DFT_Vec & y, double &dfda)
{
  double f = 0;

  x_.set(y,d_,alf+dalf);

  bool OK = true;
  do {
    try {
      f = getDF_DX();
      OK = true;
    } catch(Eta_Too_Large_Exception &e) {
      dalf /= 2;
      x_.set(y,d_,alf+dalf);
      OK = false;
    }
  } while(!OK);
  dfda = dF_.dotWith(d_);

  return f;
}



double ConjugateGradients2::step()
{

  double delta_d = d_.dotWith(d_); //dot(d_,d_);


  cout << "delta_d = " << delta_d << endl;

  double accum = 0.0;

  DFT_Vec y(x_);

  double a0 = 0.0;
  double f0 = getDF_DX();
  double df0 = dF_.dotWith(d_); //dot(dF_,d_);
  double da0 = 0;

  double da1 = -sigma_conj_grad_secant_/sqrt(delta_d); 
  double df1 = 0;
  double f1 = tryStep(a0, da1,y, df1); 
  double a1 = a0+da1;

  // we expect to move in the down-hill direction so we need f1 < f0
  if(f0 < f1) { swap(f0,f1); swap(df0,df1); swap(a0,a1);   da1 = -da1;}

  double da2 = da1; 
  double df2 = 0.0;
  double f2 = tryStep(a1,da2, y, df2); 
  double a2 = a1+da2;
  
  bool hitWall = false;

  while(f2 < f1 || (df0*df1 > 0 && df1*df2 > 0))
    {
      f0 = f1; df0 = df1; a0 = a1; da0 = da1; 
      f1 = f2; df1 = df2; a1 = a2; da1 = da2; 
      da2 *= 1.4;
      f2 = tryStep(a1,da2,y, df2);
      a2 = a1 + da2; 

      if(fabs(da2) < 1e-10) 
	{
	  hitWall = true;
	  break;
	}    
      cout << "\t\t\tbracket: " << a0 << " " << f0 << " " << df0 << " | " << a1 << " " << f1 << " " << df1 << " | " << a2 << " " << f2 << " " << df2 << endl;
    }

  if(hitWall)
    throw std::runtime_error("Hit wall while bracketing - maybe need more points on sphere?");

  cout << endl;
  cout << "\tBracketed a0 a1 a2: " << a0 << " " << a1 << " " << a2 << endl;
  cout << "\tBracketed f0 f1 f2: " << f0 << " " << f1 << " " << f2 << " : " << (f0 > f1) << " " << (f2 > f1) << endl;
  cout << "\tBracketed df0 df1 df2: " << df0 << " " << df1 << " " << df2 << endl;
  cout << endl;

  if(df1*df0 > 0 && df1*df2 > 0) 
    {
      cout.precision(12);
      cout << endl;
      cout << "\tBracketed a0 a1 a2: " << a0 << " " << a1 << " " << a2 << endl;
      cout << "\tBracketed f0 f1 f2: " << f0 << " " << f1 << " " << f2 << " : " << (f0 > f1) << " " << (f2 > f1) << endl;
      cout << "\tBracketed df0 df1 df2: " << df0 << " " << df1 << " " << df2 << endl;
      cout << endl;
	  
      throw std::runtime_error("Could not bracket");
      d_.set(r_); // = r_;
      return F_;
    }

  if(std::isnan(f0) || std::isnan(f1) || std::isnan(f2))
    {
      d_.set(r_); // = r_;
      return F_;
    }



  double xa = a0;
  double dfa = df0;
  double fa = f0;

  double xb = a1;
  double dfb = df1;
  double fb = f1;

  if(df0*df1 > 0) 
    {
      xa = a2;
      dfa = df2;
      fa = f2;      
    }

  double limit = min(min(fabs(dfa), fabs(dfb))/100, df_limit_conj_grad_); 

  for(int ix = 0; ix < jmax_conj_grad_secant_; ix++)
    {
     
      //      double x2 = xa -(ya/(yb-ya))*(xb-xa); //(xa+xb)/2;
      double x3 = (xa+xb)/2;
      //      x_ = y + x3*d_;
      x_.set(y,d_,x3);
      double f3 = getDF_DX();
      double df3 = dF_.dotWith(d_);

      if(df3*dfa > 0) { xa = x3; dfa = df3; fa = f3;}
      else           { xb = x3; dfb = df3; fb = f3;}


      cout << "\t\t\txa = " << xa << " xb = " << xb << " x3 = " << x3 << endl;
      cout << "\t\t\tfa = " << fa << " fb = " << fb << " f3 = " << f3 << endl;
      cout << "\t\t\tdfa = " << dfa << " dfb = " << dfb << " df3 = " << df3 << endl;
      cout << endl;

      if(max(fabs(dfa), fabs(dfb)) < limit) break;
    }
  
  x_.set(y,d_,Mixing_*(fa < fb ? xa : xb));

  F_ = getDF_DX();


  draw_during();  

  string s("snapshot.dat");
  density_.writeDensity(s);


  double delta_old = delta_new_;
  double delta_mid = -dF_.dotWith(r_); 

  r_.set(dF_); r_.multBy(-1); //r_ = -dF_;
  delta_new_ = r_.dotWith(r_);
  double beta = (delta_new_-delta_mid)/delta_old;

  if((step_counter_ > 0 && step_counter_%n_reset_conj_grad_ == 0) || beta <= 0)
    {
      d_.set(r_); 
      cout << "Normal reset of conjugate gradients" << endl;
    } else {d_.multBy(beta); d_.addTo(r_);} //d_ = r_+beta*d_; 


  cout << "delta_old = " << delta_old << " delta_new = " << delta_new_ << " delta_mid = " << delta_mid << " beta = " << beta << endl;

  cout << endl;

  return F_;
}

void ConjugateGradients2_Fixed::draw_after()
{
  cout << "After conj grad step " << step_counter_ << " F = " << F_ << " N = " << density_.getNumberAtoms() << " calls = " << calls_ << " delta_new_/delta_0_ = " << delta_new_/delta_0_ << " calls = " << calls_ << endl;
}

int ConjugateGradients2_Fixed::draw_during()
{
  cout << "\t1D minimization F = " << F_  << " mu_eff = " << mu_eff_ << " F-mu*N = " << F_-mu_eff_*N_fixed_target_ << " N = " << density_.getNumberAtoms() << " Ntarget = " << N_fixed_target_ << endl;
}

// The model is 
// density[i] = SMALL_VALUE + a*x[i]*x[i]
// with a chosen so that the total number is fixed:
//    N0 = (num entries)*SMALL_VALUE*dV
//    a = (N-N0)/sum(x[i]*x[i]*dV)

double ConjugateGradients2_Fixed::getDF_DX()
{
  calls_++;

  double dV = density_.dV();

  double N0 = x_.size()*SMALL_VALUE*dV;

  double sum = dV*x_.dotWith(x_); 

  double a = (N_fixed_target_-N0)/sum;

  //  vec y = sqrt(a)*x_;

  DFT_Vec y(x_);
  y.multBy(sqrt(a));

  density_.set_density_from_amplitude(y);

  double F = 0;

  try {
    F = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF_,false);   
    cout.precision(12);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  mu_eff_ = 0.0;
  for(long i=0;i<dF_.size();i++)
    mu_eff_ += dF_.get(i)*a*x_.get(i)*x_.get(i);
  mu_eff_ *= dV/(N_fixed_target_-N0);

  //  dF_.addTo(-mu_eff_);


  double ff = 0;
  double im = 0;
  double dm = 0;
  for(int i=0;i<dF_.size();i++)
    {
      double force = dF_.get(i)-mu_eff_;

      double density = density_.getDensity(i);
      double density_predicted = exp(log(density)-force/dV);
      if(fabs(density-density_predicted) > ff) 
	{ ff = fabs(density-density_predicted); im = i;}
      dF_.set(i,2*a*force*x_.get(i));
    }

  err_ = ff;

  //  dF_.Schur(x_,dF_);
  //  dF_.multBy(2*a);


  mu1 = mu_eff_;

  return F;
}

/*
void ConjugateGradients2::check(bool write)
{
  double F = 0;
  vec dF(density_.Ntot());
  F = getDF_DX(); //dft_.calculateFreeEnergyAndDerivatives(density_,mu_, dF,false);   
  dF = dF_;

  double N = density_.getNumberAtoms()/density_.dV();


  double eps = 1e-5;
  vec dF1(density_.Ntot());

  double biggest = 0;
  double dmin = 100;
  int ib,jb,kb;

  ofstream db("db.dat");
  cout << "mu_ = " << mu_ << endl;
  for(int i=37;i<density_.Nx();i++)
    for(int j=63;j<density_.Ny();j++)
      for(int k=63;k<density_.Nz();k++)
	{
	  long pos = density_.pos(i,j,k);
	  double d = x_(pos); 

	  if(d < 0.5) continue;
	  cout << i << " " << j << " " << k << endl;

	  if(d < dmin) {cout << "\t\t\t\t\t\t\t\t\t\t\ti = " << i << " j = " << j << " k = " << k << " dmin = " << d << endl; dmin = d;}

	  double dp = d*(1+eps);
	  x_(pos) = dp;
	  density_.set_density_from_amplitude(x_); 
	  double Fp = getDF_DX();
	  double Np = density_.getNumberAtoms();

	  double dm = d*(1-eps);
	  x_(pos) = dm;

	  double Fm = getDF_DX();
	  double Nm = density_.getNumberAtoms();


	  x_(pos) = d;

	  double est = (Fp-Fm)/(2*eps*d);

	  if(fabs(est-dF[pos]) > 0.1*fabs(est) || (i == 37 && j == 63 && k == 63)) //biggest)
	    {

	      double dsum = 0;
	      for(long pos = 0; pos < density_.Ntot(); pos++)
		dsum += density_.getDensity(pos);
	      

	      biggest = fabs(est-dF[pos]);
	      ib = i;
	      jb = j;
	      kb = k;

	      cout.precision(12);

	      cout << i << " " << j << " " << k << " " << pos << " " << d << " " << dF[pos] << " " << (Fp-Fm)/(2*eps*d) << " " << Np << " " << Nm << " " << (Np-Nm)/(2*eps*d) << endl;
	      cout << "\t\t" << dft_.getEta(pos,0) << " " << dsum << " " << N << " " << Np << " " << Nm << " " << Fp << " " << Fm << " " << endl;
	    }
	  if(write) db << i << " " << j << " " << k << " " << d << " " << dF[pos] << " " << (Fp-Fm)/(2*eps*d) << endl;
	}
  db.close();


}
*/

/*
int nloptCount = 0;
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
  ++nloptCount;

  nlOptMinimizer &opt = *((nlOptMinimizer*) my_func_data);

  return opt.objectiveFunction(n,x,grad);
}

double nlOptMinimizer::objectiveFunction(unsigned n, const double *x, double *grad)
{  
  cout << "n = " << n << endl;

  double dV = density_.dV();
  double N0 = n*SMALL_VALUE*dV;

  DFT_Vec y;   y.set(x,n);
  
  double sum2 = y.dotWith(y); 
  double a = (N_fixed_target_-N0)/(sum2*dV);
  y.multBy(sqrt(a));
  density_.set_density_from_amplitude(y);
  
  double F = 0;
  DFT_Vec dF(n);


  try {
    F = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF,false);   
    cout.precision(12);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  cout << "F = " << F;
  
  if(grad)
    {
      double mu_eff = 0.0;
  	
      for(long i=0;i<dF.size();i++)
	mu_eff += dF.get(i)*x[i]*x[i];
      mu_eff /= sum2;

      double mm = 0;
      
      for(int i=0;i<dF.size();i++)
	{
	  grad[i] = 2*a*x[i]*(dF.get(i)-mu_eff)/dV; // because dF carries a spare factor of dV.
	  mm = max(mm,grad[i]);
	}
      cout << " dfmax = " << mm;      
    }
  cout << endl;
  return F;
}


void nlOptMinimizer::run(string& logfile, long maxSteps)
{
  nlopt::opt opt(nlopt::LN_BOBYQA, density_.Ntot());
  //    nlopt::opt opt(nlopt::LD_LBFGS, density_.Ntot());
  //  opt.set_lower_bounds(SMALL_VALUE);
  opt.set_min_objective(myfunc, (void*) this);
  opt.set_xtol_rel(1e-4);

  std::vector<double> x(density_.Ntot());
  for(long i=0;i<density_.Ntot();i++)
    x[i] = sqrt(max(0.0, density_.getDensity(i) - SMALL_VALUE));

  double minf;
  nlopt::result result = opt.optimize(x, minf);  
}
*/


void fireMinimizer::initialize()
{
  Minimizer_Fixed_N::initialize();
  
  it_ = 0;
  cut_ = 0;


  alpha_ = alpha_start_;
  
  F_ = getDF_DX();
  cout << "Here: F_ = " << F_ << endl;
  v_.zeros(v_.size());
}


void fireMinimizer::verlet()
{
  DFT_Vec y(x_);
  //  DFT_Vec yv(v_);
  
  try{
    // this gives x_(t+dt)
    x_.Increment_And_Scale(v_,dt_);
    x_.Increment_And_Scale(dF_,-0.5*dt_*dt_); // dF=dV/dx does not have minus sign ... 
    // now do half the v update
    v_.Increment_And_Scale(dF_, -0.5*dt_);
    // then get new forces
    F_ = getDF_DX();
    // and finish velocity update
    v_.Increment_And_Scale(dF_, -0.5*dt_);
  } catch (Eta_Too_Large_Exception &e) {
    x_.set(y);
    //    v_.set(yv);
    v_.zeros(v_.size());
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
  
  string s("snapshot.dat");
  density_.writeDensity(s);

  
  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  double p = -v_.dotWith(dF_);

  v_.multBy(1-alpha_);
  v_.Increment_And_Scale(dF_, -alpha_*v_.euclidean_norm()/dF_.euclidean_norm());

  if(p < 0)
    {
      cout << "\tp < 0 : reset" << endl;
      v_.zeros(v_.size());
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
  cout << "After FIRE step " << step_counter_ << " F = " << F_ << " N = " << density_.getNumberAtoms() << " calls = " << calls_ << " dt_ = " << dt_ << " alpha = " << alpha_ << endl;
}

int fireMinimizer::draw_during()
{
  cout << "\t1D minimization F = " << F_  << " mu_eff = " << mu_eff_ << " F-mu*N = " << F_-mu_eff_*N_fixed_target_ << " N = " << density_.getNumberAtoms() << " Ntarget = " << N_fixed_target_ << endl;
}


//Here
void fireMinimizer_Mu::initialize()
{
  Minimizer::initialize();
  
  it_ = 0;
  cut_ = 0;

  alpha_ = alpha_start_;
  
  F_ = getDF_DX();
  cout << "Here: F_ = " << F_ << endl;
  v_.zeros(v_.size());
}


void fireMinimizer_Mu::verlet()
{
  DFT_Vec y(x_);
  //  DFT_Vec yv(v_);
  
  try{
    // this gives x_(t+dt)
    x_.Increment_And_Scale(v_,dt_);
    x_.Increment_And_Scale(dF_,-0.5*dt_*dt_); // dF=dV/dx does not have minus sign ... 
    // now do half the v update
    v_.Increment_And_Scale(dF_, -0.5*dt_);
    // then get new forces
    F_ = getDF_DX();
    // and finish velocity update
    v_.Increment_And_Scale(dF_, -0.5*dt_);
  } catch (Eta_Too_Large_Exception &e) {
    x_.set(y);
    //    v_.set(yv);
    v_.zeros(v_.size());
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
  
  string s("snapshot.dat");
  density_.writeDensity(s);

  
  // dF does not include the minus so we have to put it in by hand everywhere from here down:
  double p = -v_.dotWith(dF_);

  v_.multBy(1-alpha_);
  v_.Increment_And_Scale(dF_, -alpha_*v_.euclidean_norm()/dF_.euclidean_norm());

  if(p < 0)
    {
      v_.zeros(v_.size());
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
  cout << "After FIRE step " << step_counter_ << " F = " << F_ << " N = " << density_.getNumberAtoms() << " calls = " << calls_ << " dt_ = " << dt_ << " alpha = " << alpha_ << endl;
}

int fireMinimizer_Mu::draw_during()
{
  cout << "\t1D minimization F-mu*N = " << F_  << " N = " << density_.getNumberAtoms() << endl;
}
