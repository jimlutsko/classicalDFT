#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>

#include "Spliner.h"

using namespace std;


Spliner::Spliner(double const * const x, double const * const y, int N, double y0, double yn)
  : N_(N)
{
  // copy data
  x_ = new double[N_];
  memcpy(x_,x,N_*sizeof(double));
  y_ = new double[N_];
  memcpy(y_,y,N_*sizeof(double));

  // Set up the cubic spline. This is basically the NR spline(...) routine.

  y2_ = new double[N_];

  double *u = new double[N_-1];
  if (y0 > 0.99e30)
    y2_[0]=u[0]=0.0;
  else {
    y2_[0] = -0.5;
    u[0]=(3.0/(x_[1]-x_[0]))*((y_[1]-y_[0])/(x_[1]-x_[0])-y0);
  }
  for (int i=1;i<N_-1;i++) 
    {
      double sig=(x_[i]-x_[i-1])/(x_[i+1]-x_[i-1]);
      double p=sig*y2_[i-1]+2.0;
      y2_[i]=(sig-1.0)/p;
      u[i]=(y_[i+1]-y_[i])/(x_[i+1]-x_[i]) - (y_[i]-y_[i-1])/(x_[i]-x_[i-1]);
      u[i]=(6.0*u[i]/(x_[i+1]-x_[i-1])-sig*u[i-1])/p;
    }
  double qn,un;
  if (yn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x_[N_-1]-x_[N_-2]))*(yn-(y_[N_-1]-y_[N_-2])/(x_[N_-1]-x_[N_-2]));
  }
  y2_[N_-1]=(un-qn*u[N_-2])/(qn*y2_[N_-2]+1.0);
  for (int k=N_-2;k>=0;k--)
    y2_[k]=y2_[k]*y2_[k+1]+u[k];
  delete[] u;
}

Spliner::Spliner(Spliner &s)
{
  N_ = s.N_;
  x_ = new double[N_];
  memcpy(x_,s.x_,N_*sizeof(double));
  y_ = new double[N_];
  memcpy(y_,s.y_,N_*sizeof(double));
  y2_ = new double[N_];
  memcpy(y2_,s.y2_,N_*sizeof(double));
}

Spliner::~Spliner()
{
    if(x_ != NULL) delete[] x_; x_ = NULL;
    if(y_ != NULL) delete[] y_; y_ = NULL;
    if(y2_ != NULL) delete[] y2_; y2_ = NULL;
}


double Spliner::integrate(double x1, double x2)
{
  if(x2 < x1) throw std::runtime_error("x2 < x1 in spliner.integrate");



  int klo1, khi1;
  double a1,b1,h1;
  locate(x1,a1,b1,h1, klo1, khi1);

  int klo2, khi2;
  double a2,b2,h2;
  locate(x2,a2,b2,h2, klo2, khi2);

  //nb a = x-xlo  and  b = xhi-x

  double ret = 0.0;

  if(klo1 == klo2) // cubic integral from x1 to x2
  {
      ret = -0.5*h1*(a2*a2-a1*a1)*y_[klo1]+0.5*h1*(b2*b2-b1*b1)*y_[khi1]+(-(0.25*h1*(a2*a2*a2*a2-a1*a1*a1*a1)-0.5*h1*(a2*a2-a1*a1))*y2_[klo1]
									  +(0.25*h1*(b2*b2*b2*b2-b1*b1*b1*b1)-0.5*h1*(b2*b2-b1*b1))*y2_[khi1])*(h1*h1)/6.0;
  } else {
      // contribution of x1 to x[khi1]
      ret = -0.5*h1*(0-a1*a1)*y_[klo1]+0.5*h1*(1-b1*b1)*y_[khi1]+(-(0.25*h1*(0-a1*a1*a1*a1)-0.5*h1*(0-a1*a1))*y2_[klo1]
								  +(0.25*h1*(1-b1*b1*b1*b1)-0.5*h1*(1-b1*b1))*y2_[khi1])*(h1*h1)/6.0;
      // contribution of intervening intervals
      for(int k = khi1; k < klo2; k++)
      {
	  double h = x_[k+1]-x_[k];
	  ret += 0.5*h*(y_[k]+y_[k+1])-(y2_[k]+y2_[k+1])*(h*h*h)/24.0;
      }
      // contribution of x[klo2] to x2
      ret += -0.5*h2*(a2*a2-1)*y_[klo2]+0.5*h2*(b2*b2-0)*y_[khi2]+(-(0.25*h2*(a2*a2*a2*a2-1)-0.5*h2*(a2*a2-1))*y2_[klo2]
								  +(0.25*h2*(b2*b2*b2*b2-0)-0.5*h2*(b2*b2-0))*y2_[khi2])*(h2*h2)/6.0;

  }
  return ret;
}



double Spliner::f(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  double yz =a*y_[klo]+b*y_[khi]+((a*a*a-a)*y2_[klo]
			 +(b*b*b-b)*y2_[khi])*(h*h)/6.0;
  return yz;
}

double Spliner::dfdx(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  double yz =-y_[klo]+y_[khi]+(-(3*a*a-1)*y2_[klo]
		      +(3*b*b-1)*y2_[khi])*(h*h)/6.0;
  yz /= h;
  return yz;
}

double Spliner::d2fdx2(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  return a*y2_[klo]+b*y2_[khi];
}

double Spliner::d3fdx3(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  return (y2_[khi]-y2_[klo])/h;
}

void Spliner::locate(double z,  double &a, double &b, double &h, int &klo, int &khi) const
{
  klo=0;
  khi=N_-1;
  while (khi-klo > 1) 
    {
      int k=(khi+klo) >> 1;
      if (x_[k] > z) khi=k;
      else klo=k;
    }
  h=x_[khi]-x_[klo];
  if (h == 0.0) 
    throw std::runtime_error("Bad xa input to Spliner::locate");
  a=(x_[khi]-z)/h;
  b=(z-x_[klo])/h;
}


Spliner2::Spliner2(double const * const x1, 
		   double const * const * const x2, 
		   double const * const * const y, 
		   int N1, int N2, double bc)
  : N1_(N1), N2_(N2), bc_(bc)
{
  construct(x1,x2,y);
}

void Spliner2::construct(double const * const x1, 
		    double const * const * const x2, 
		    double const * const * const y)
{

  x1_ = new double[N1_];
  memcpy(x1_,x1,N1_*sizeof(double));

  s1_ = new Spliner*[N1_];
  for (int j=0; j<N1_;j++) 
    s1_[j] = new Spliner(x2[j], y[j], N2_, 1e30, 1e30);

  // working space
  int NMAX = (N1_ > N2_ ? N1_ : N2_);
  ytmp_  = new double[NMAX];
  yytmp_ = new double[NMAX];
}

Spliner2::Spliner2(const string &file, double bc)
{
  ifstream in(file.c_str());

  if(in.bad())
      throw std::runtime_error("Matrix file  cannot be opened ... aborting!");

  int Nr, Nx;

  Nr = Nx = 0;

  char buf[256];
  int newRow = 1;
  while(!in.eof())
    {
      in.getline(buf,256);
      if(buf[0])
	{
	  if(buf[0] == '#') continue;	 
	  if(newRow){ newRow = 0; Nx = 1;}
	  else Nx++;
	} else {Nr++; newRow = 1;}
    }
  in.close();


  double *r = new double[Nr];  
  double **x = new double*[Nr];
  x[0] = new double[(Nr+1)*Nx];
  double **phi = new double*[Nr];
  phi[0] = new double[(Nr+1)*(Nx)];

  for(int i=1;i<Nr;i++)
    {
      x[i] = x[i-1]+Nx;
      phi[i] = phi[i-1]+Nx;
    }

  ifstream in1(file.c_str());
  int _nr, _nx;
  _nr = _nx = 0;
  newRow = 1;
  while(!in1.eof())
    {
      in1.getline(buf,256);
      if(buf[0])
	{
	  if(buf[0] == '#') 
	    continue;	 

	  char *pch = strtok (buf,"\t");
	  r[_nr] = atof(pch); 
	  pch = strtok(NULL, "\t");

	  if(newRow) {_nx = 0; newRow = 0;}
	  
	  x[_nr][_nx] = atof(pch); 
	  pch = strtok(NULL, "\t");

	  phi[_nr][_nx] = atof(pch); 

	  _nx++;

	} else {_nr++; newRow = 1;}
    }
  in1.close();  

  N1_ = Nr;
  N2_ = Nx;
  bc_ = bc;

  construct(r,x,phi);
  delete r;
  delete x;
  delete phi;
}




Spliner2::~Spliner2()
{
  delete[] x1_;
  //  delete[] (*x2_);
  //  delete[] x2_;
  //  delete[] (*y_);
  //  delete[] y_;

  if(s1_)
      for(int i=0;i<N1_; i++) 
	  delete (s1_[i]);

  delete s1_; s1_ = NULL;

  delete ytmp_; ytmp_ = NULL;
  delete yytmp_; yytmp_ = NULL;
}


double Spliner2::getMaxDensity()
{
  throw runtime_error("Spliner2::getMaxDensity() no longer implemented");
  /*
    double r = 0;
    for(int j=0;j<N1_;j++)
    if(x2_[j][N2_-1] > r) r = x2_[j][N2_-1];
    return r;
  */
}


double Spliner2::f(double z1, double z2) const
{
  for (int j=0;j<N1_;j++) 
    yytmp_[j] = s1_[j]->f(z2);

  Spliner s(x1_, yytmp_, N1_, bc_, 1e30);
  return s.f(z1);
}

double Spliner2::dfdx1(double z1, double z2) const
{
  for (int j=0;j<N1_;j++) 
    yytmp_[j] = s1_[j]->f(z2);

  Spliner s(x1_, yytmp_, N1_, bc_, 1e30);
  return s.dfdx(z1);
}

double Spliner2::dfdx2(double z1, double z2) const
{
  for (int j=0;j<N1_;j++) 
    yytmp_[j] = s1_[j]->dfdx(z2);

  Spliner s(x1_, yytmp_, N1_, bc_, 1e30);
  return s.f(z1);
}

double Spliner2::d2fdx1dx1(double z1, double z2) const
{
  for (int j=0;j<N1_;j++) 
    yytmp_[j] = s1_[j]->f(z2);

  Spliner s(x1_, yytmp_, N1_, bc_, 1e30);
  return s.d2fdx2(z1);
}

double Spliner2::d2fdx2dx2(double z1, double z2) const
{
  for (int j=0;j<N1_;j++) 
    yytmp_[j] = s1_[j]->d2fdx2(z2); 

  Spliner s(x1_, yytmp_, N1_, bc_, 1e30);
  return s.f(z1);
}

double Spliner2::d2fdx1dx2(double z1, double z2) const
{
  for (int j=0;j<N1_;j++) 
    yytmp_[j] = s1_[j]->dfdx(z2);

  Spliner s(x1_, yytmp_, N1_, bc_, 1e30);
  return s.dfdx(z1);
}

///////////////// Vector based spliner


SplinerVec::SplinerVec(const vector<double> &x, const vector<double> &y, double y0, double yn)
{initialize(x,y,y0,yn);}

void SplinerVec::initialize(const vector<double> &x, const vector<double> &y, double y0, double yn)
{
    // copy ata
    x_.clear(); // because I am paranoid
    y_.clear();
    y2_.clear();

    x_ = x;
    y_ = y;

    // Set up the cubic spline. This is basically the NR spline(...) routine.

    vector<double> u;
    if (y0 > 0.99e30){
	y2_.push_back(0.0); 
	u.push_back(0.0);
    } else {
	y2_.push_back(-0.5); 
	u.push_back((3.0/(x_[1]-x_[0]))*((y_[1]-y_[0])/(x_[1]-x_[0])-y0));
    }

    double uprev = u[0];

    for (unsigned int i=1;i<x_.size()-1;i++) 
    {
	double sig=(x_[i]-x_[i-1])/(x_[i+1]-x_[i-1]);
	double p=sig*y2_[i-1]+2.0;

	double yy = (sig-1.0)/p;
	double uu = (y_[i+1]-y_[i])/(x_[i+1]-x_[i]) - (y_[i]-y_[i-1])/(x_[i]-x_[i-1]);
	uu=(6.0*uu/(x_[i+1]-x_[i-1])-sig*uprev)/p;
	
	y2_.push_back(yy);
	u.push_back(uu);
	uprev = uu;
    }
    double qn,un;
    int N = y_.size();
    if (yn > 0.99e30)
	qn=un=0.0;
    else {
	qn=0.5;
	un=(3.0/(x_[N-1]-x_[N-2]))*(yn-(y_[N-1]-y_[N-2])/(x_[N-1]-x_[N-2]));
    }
    y2_.push_back((un-qn*u[N-2])/(qn*y2_[N-2]+1.0));
    for (int k=N-2;k>=0;k--)
	y2_[k]=y2_[k]*y2_[k+1]+u[k];
}

SplinerVec::SplinerVec(SplinerVec const &s)
    : x_(s.x_), y_(s.y_), y2_(s.y2_) {}

SplinerVec::~SplinerVec(){}



double SplinerVec::integrate(double x1, double x2)
{
  if(x2 < x1) throw std::runtime_error("x2 < x1 in spliner.integrate");

  int klo1, khi1;
  double a1,b1,h1;
  locate(x1,a1,b1,h1, klo1, khi1);

  int klo2, khi2;
  double a2,b2,h2;
  locate(x2,a2,b2,h2, klo2, khi2);

  //nb a = x-xlo  and  b = xhi-x

  double ret = 0.0;

  if(klo1 == klo2) // cubic integral from x1 to x2
  {
      ret = -0.5*h1*(a2*a2-a1*a1)*y_[klo1]+0.5*h1*(b2*b2-b1*b1)*y_[khi1]
	  +(-(0.25*h1*(a2*a2*a2*a2-a1*a1*a1*a1)-0.5*h1*(a2*a2-a1*a1))*y2_[klo1]
									  
	    +(0.25*h1*(b2*b2*b2*b2-b1*b1*b1*b1)-0.5*h1*(b2*b2-b1*b1))*y2_[khi1])*(h1*h1)/6.0;
  } else {
      // contribution of x1 to x[khi1]
      ret = -0.5*h1*(0-a1*a1)*y_[klo1]+0.5*h1*(1-b1*b1)*y_[khi1]
	  +(-(0.25*h1*(0-a1*a1*a1*a1)-0.5*h1*(0-a1*a1))*y2_[klo1]
	    +(0.25*h1*(1-b1*b1*b1*b1)-0.5*h1*(1-b1*b1))*y2_[khi1])*(h1*h1)/6.0;
      // contribution of intervening intervals
      for(int k = khi1; k < klo2; k++)
      {
	  double h = x_[k+1]-x_[k];
	  ret += 0.5*h*(y_[k]+y_[k+1])-(y2_[k]+y2_[k+1])*(h*h*h)/24.0;
      }
      // contribution of x[klo2] to x2
      ret += -0.5*h2*(a2*a2-1)*y_[klo2]+0.5*h2*(b2*b2-0)*y_[khi2]
	  +(-(0.25*h2*(a2*a2*a2*a2-1)-0.5*h2*(a2*a2-1))*y2_[klo2]
	    +(0.25*h2*(b2*b2*b2*b2-0)-0.5*h2*(b2*b2-0))*y2_[khi2])*(h2*h2)/6.0;

  }
  return ret;
}



double SplinerVec::f(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  double yz =a*y_[klo]+b*y_[khi]+((a*a*a-a)*y2_[klo]
			 +(b*b*b-b)*y2_[khi])*(h*h)/6.0;
  return yz;
}

double SplinerVec::dfdx(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  double yz =-y_[klo]+y_[khi]+(-(3*a*a-1)*y2_[klo]
		      +(3*b*b-1)*y2_[khi])*(h*h)/6.0;
  yz /= h;
  return yz;
}

double SplinerVec::d2fdx2(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  return a*y2_[klo]+b*y2_[khi];
}

double SplinerVec::d3fdx3(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  return (y2_[khi]-y2_[klo])/h;
}

void SplinerVec::locate(double z,  double &a, double &b, double &h, int &klo, int &khi) const
{
  klo=0;
  khi=x_.size()-1;
  while (khi-klo > 1) 
    {
      int k=(khi+klo) >> 1;
      if (x_[k] > z) khi=k;
      else klo=k;
    }
  h=x_[khi]-x_[klo];
  if (h == 0.0) 
    throw std::runtime_error("Bad xa input to SplinerVec::locate");
  a=(x_[khi]-z)/h;
  b=(z-x_[klo])/h;
}



///////////////// Pair based spliner


void SplinerPair::initialize(vector< pair<double,double> > const &data, double y0, double yn)
{
    // copy data
    data_.clear(); // because I am paranoid
    y2_.clear();

    data_ = data;
    // Set up the cubic spline. This is basically the NR spline(...) routine.

    vector<double> u;

    if (y0 > 0.99e30){
	y2_.push_back(0.0); 
	u.push_back(0.0);
    } else {
	y2_.push_back(-0.5); 
	u.push_back((3.0/(x(1)-x(0)))*((y(1)-y(0))/(x(1)-x(0))-y0));
    }

    double uprev = u[0];

    for (unsigned int i=1;i<data_.size()-1;i++) 
    {
	double sig=(x(i)-x(i-1))/(x(i+1)-x(i-1));
	double p=sig*y2_[i-1]+2.0;

	double yy = (sig-1.0)/p;
	double uu = (y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1));
	uu=(6.0*uu/(x(i+1)-x(i-1))-sig*uprev)/p;
	
	y2_.push_back(yy);
	u.push_back(uu);
	uprev = uu;
    }
    double qn,un;
    int N = data_.size();
    if (yn > 0.99e30)
	qn=un=0.0;
    else {
	qn=0.5;
	un=(3.0/(x(N-1)-x(N-2)))*(yn-(y(N-1)-y(N-2))/(x(N-1)-x(N-2)));
    }
    y2_.push_back((un-qn*u[N-2])/(qn*y2_[N-2]+1.0));
    for (int k=N-2;k>=0;k--)
	y2_[k]=y2_[k]*y2_[k+1]+u[k];
}

SplinerPair::SplinerPair(SplinerPair const &s)
    : data_(s.data_), y2_(s.y2_) {}

SplinerPair::~SplinerPair(){}



double SplinerPair::integrate(double x1, double x2)
{
  if(x2 < x1) throw std::runtime_error("x2 < x1 in spliner.integrate");

  int klo1, khi1;
  double a1,b1,h1;
  locate(x1,a1,b1,h1, klo1, khi1);

  int klo2, khi2;
  double a2,b2,h2;
  locate(x2,a2,b2,h2, klo2, khi2);

  //nb a = x-xlo  and  b = xhi-x

  double ret = 0.0;

  if(klo1 == klo2) // cubic integral from x1 to x2
  {
      ret = -0.5*h1*(a2*a2-a1*a1)*y(klo1)+0.5*h1*(b2*b2-b1*b1)*y(khi1)
	  +(-(0.25*h1*(a2*a2*a2*a2-a1*a1*a1*a1)-0.5*h1*(a2*a2-a1*a1))*y2_[klo1]
									  
	    +(0.25*h1*(b2*b2*b2*b2-b1*b1*b1*b1)-0.5*h1*(b2*b2-b1*b1))*y2_[khi1])*(h1*h1)/6.0;
  } else {
      // contribution of x1 to x[khi1]
      ret = -0.5*h1*(0-a1*a1)*y(klo1)+0.5*h1*(1-b1*b1)*y(khi1)
	  +(-(0.25*h1*(0-a1*a1*a1*a1)-0.5*h1*(0-a1*a1))*y2_[klo1]
	    +(0.25*h1*(1-b1*b1*b1*b1)-0.5*h1*(1-b1*b1))*y2_[khi1])*(h1*h1)/6.0;
      // contribution of intervening intervals
      for(int k = khi1; k < klo2; k++)
      {
	  double h = x(k+1)-x(k);
	  ret += 0.5*h*(y(k)+y(k+1))-(y2_[k]+y2_[k+1])*(h*h*h)/24.0;
      }
      // contribution of x[klo2] to x2
      ret += -0.5*h2*(a2*a2-1)*y(klo2)+0.5*h2*(b2*b2-0)*y(khi2)
	  +(-(0.25*h2*(a2*a2*a2*a2-1)-0.5*h2*(a2*a2-1))*y2_[klo2]
	    +(0.25*h2*(b2*b2*b2*b2-0)-0.5*h2*(b2*b2-0))*y2_[khi2])*(h2*h2)/6.0;

  }
  return ret;
}



double SplinerPair::f(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  double yz =a*y(klo)+b*y(khi)+((a*a*a-a)*y2_[klo]
			 +(b*b*b-b)*y2_[khi])*(h*h)/6.0;
  return yz;
}

double SplinerPair::dfdx(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  double yz =-y(klo)+y(khi)+(-(3*a*a-1)*y2_[klo]
		      +(3*b*b-1)*y2_[khi])*(h*h)/6.0;
  yz /= h;
  return yz;
}

double SplinerPair::d2fdx2(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  return a*y2_[klo]+b*y2_[khi];
}

double SplinerPair::d3fdx3(double z) const 
{
  int klo, khi;
  double a,b,h;
  locate(z,a,b,h, klo, khi);
  return (y2_[khi]-y2_[klo])/h;
}

void SplinerPair::locate(double z,  double &a, double &b, double &h, int &klo, int &khi) const
{
  klo=0;
  khi=data_.size()-1;
  while (khi-klo > 1) 
    {
      int k=(khi+klo) >> 1;
      if (x(k) > z) khi=k;
      else klo=k;
    }
  h=x(khi)-x(klo);
  if (h == 0.0) 
    throw std::runtime_error("Bad xa input to SplinerPair::locate");
  a=(x(khi)-z)/h;
  b=(z-x(klo))/h;
}

///////////////////////////////////////////////
// General 2-D spliner

Spliner2General::Spliner2General(vector< pair< double, vector< pair<double,double> > > > const &data,double bc, double bc20,double bc21)
{
    initialize(data,bc,bc20,bc21);
}

void Spliner2General::initialize(vector< pair< double, vector< pair<double,double> > > >  const &data,double bc, double bc20, double bc21)
{
    data_ = data; 
    bc_ = bc;

    for (unsigned j=0; j<data.size();j++) 
    {

	vector< pair<double,double> > v = data_[j].second;
	pair<double,double> p1 = v[v.size()-2]; 
	pair<double,double> p2 = v[v.size()-1]; 
	double slope = ((p2.second-p1.second)/(p2.first-p1.first))*((p2.second-p1.second)/(p2.first-p1.first));



	SplinerPair s(data_[j].second,bc20,slope);
	s1_.push_back(s);
    }
}


double Spliner2General::getMaxDensity()
{
  throw runtime_error("Spliner2General::getMaxDensity() no longer implemented");
}


double Spliner2General::f(double z1, double z2) const
{
    vector< pair <double, double> > yytmp;

    for (unsigned j=0;j<data_.size();j++)
    {
	pair<double, double>  p;
	p.first  = data_[j].first;
	p.second = s1_[j].f(z2);
	yytmp.push_back(p);
    }
    SplinerPair s(yytmp, bc_, 1e30);
    return s.f(z1);
}


double Spliner2General::dfdx1(double z1, double z2) const
{
    vector< pair <double, double> > yytmp;
    
    for (unsigned j=0;j<data_.size();j++)
    {
	pair<double, double>  p;
	p.first  = data_[j].first;
	p.second = s1_[j].f(z2);
	yytmp.push_back(p);
    }
    SplinerPair s(yytmp, bc_, 1e30);

    return s.dfdx(z1);
}

double Spliner2General::dfdx2(double z1, double z2) const
{
    vector< pair <double, double> > yytmp;
    
    for (unsigned j=0;j<data_.size();j++)
    {
	pair<double, double>  p;
	p.first  = data_[j].first;
	p.second = s1_[j].dfdx(z2);
	yytmp.push_back(p);
    }
    SplinerPair s(yytmp, bc_, 1e30);

    return s.f(z1);
}

double Spliner2General::d2fdx1dx1(double z1, double z2) const
{
    vector< pair <double, double> > yytmp;
    
    for (unsigned j=0;j<data_.size();j++)
    {
	pair<double, double>  p;
	p.first  = data_[j].first;
	p.second = s1_[j].f(z2);
	yytmp.push_back(p);
    }
    SplinerPair s(yytmp, bc_, 1e30);
    return s.d2fdx2(z1);
}

double Spliner2General::d2fdx2dx2(double z1, double z2) const
{
    vector< pair <double, double> > yytmp;
    
    for (unsigned j=0;j<data_.size();j++)
    {
	pair<double, double>  p;
	p.first  = data_[j].first;
	p.second = s1_[j].d2fdx2(z2);
	yytmp.push_back(p);
    }
    SplinerPair s(yytmp, bc_, 1e30);

    return s.f(z1);
}

double Spliner2General::d2fdx1dx2(double z1, double z2) const
{
    vector< pair <double, double> > yytmp;
    
    for (unsigned j=0;j<data_.size();j++)
    {
	pair<double, double>  p;
	p.first  = data_[j].first;
	p.second = s1_[j].dfdx(z2);
	yytmp.push_back(p);
    }
    SplinerPair s(yytmp, bc_, 1e30);
    return s.dfdx(z1);
}

