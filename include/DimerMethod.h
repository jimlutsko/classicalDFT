#ifndef __LUTSKO__DIMERMETHOD__
#define __LUTSKO__DIMERMETHOD__



#include "DFT.h"

class DimerMethod
{
 public:

 DimerMethod(Density &density, DFT &dft, DFT_Vec &nvec, double eps, bool verbose = false)
   : density_(density), dft_(dft), R0_(density.size()), R1_(density.size()), R2_(density.size()), N_(nvec), T_(density.size()), dR_(eps)
  {      

    for(long i=0;i<density_.Ntot();i++)
      R0_.set(i,sqrt(max(0.0, density_.getDensity(i) - SMALL_VALUE)));

    F0_.zeros(R0_.size());
    F1_.zeros(R0_.size());
    F2_.zeros(R0_.size());

    V_.zeros(R0_.size());
      
    verbose_ = verbose;

    dV_ = density.dV();
    Nfixed_ = density.getNumberAtoms();

    R1_.set(R0_);  R1_.Increment_And_Scale(N_,dR_);
    R2_.set(R0_);  R2_.Increment_And_Scale(N_,-dR_);
  }
  
  ~DimerMethod(){}

  double getR0(int i) { return (R1_.get(i)+R2_.get(i))/2;}
  double getR1(int i) { return R1_.get(i);}
  double getR2(int i) { return R2_.get(i);}


  bool findTransitionState(double dt = 0.0001, long itmax = 10000)
  {
    // Initialization
    EnergyAndForce(R0_,E0_,F0_);    
    
    for(int i=0;i<itmax;i++)
      {
	if(verbose_) cout << "step = " << i << endl;
	rotateDimer(false);
	if(translateDimer(dt)) {cout << "step = " << i << endl; return true;}
	if(verbose_) cout << endl;           
      }
    return false;
  }

  void rotateDimer(bool doCheck = false)
  {
    // Get the initial energies and forces
    EnergyAndForce(R1_,E1_,F1_);
    /*
    // Determine second vector that defines the rotational plane
    // This is just the component of F1-F2 perpindicular to N    
    // F2 = 2F0-F1
    // T = F1-F2 = 2(F1-F0)

    T_.set(F1_);
    T_.DecrementBy(F0_);

    double F10_dot_N = T_.dotWith(N_);
    T_.Increment_And_Scale(N_,-F10_dot_N);
    double T_norm = T_.euclidean_norm(); // we need this below
    T_.normalise();
    
    cout << "Finished constructing Theta for phi = 0" << endl;

    double F12_dot_N = 2*F10_dot_N;
    
    // E2
    //    E2_ = 2*E0_-0.5*dR_*F12_dot_N-E1_; 

    double E0 = 2*E0_-0.5*dR_*F12_dot_N;
    */

    EnergyAndForce(R2_,E2_,F2_);
    T_.set(F1_);
    T_.DecrementBy(F2_);

    double F12_dot_N = T_.dotWith(N_);
    T_.Increment_And_Scale(N_,-F12_dot_N);
    double T_norm = T_.euclidean_norm()/2; // we need this below
    T_.normalise();
    double E0 = E1_+E2_;
    
    double F0 = 2*T_norm;

    // rotation angle
    double phi1 = 0.1; //0.5*atan(-dC_dPhi_0/(2*C0));

    DFT_Vec N_rotated(N_);
    N_rotated.multBy(cos(phi1));
    N_rotated.Increment_And_Scale(T_,sin(phi1));
    
    DFT_Vec R1_rotated(R0_);
    R1_rotated.Increment_And_Scale(N_rotated,dR_);

    double E1_rotated;
    DFT_Vec F1_rotated; F1_rotated.zeros(F1_.size());
    EnergyAndForce(R1_rotated,E1_rotated,F1_rotated);

    //    double F12_dot_N_rotated = 2*(F1_rotated.dotWith(N_rotated)-F0_.dotWith(N_rotated));

    // Rotated Curvature and gradient
    //    double Ephi = 2*E0_-0.5*dR_*F12_dot_N_rotated;

    DFT_Vec R2_rotated(R0_);
    R2_rotated.Increment_And_Scale(N_rotated,-dR_);

    double E2_rotated;
    DFT_Vec F2_rotated; F2_rotated.zeros(F2_.size());
    EnergyAndForce(R2_rotated,E2_rotated,F2_rotated);
    double Ephi = E1_rotated+E2_rotated;

    //////////// Fourier coefficients
    double C = -0.5*dR_*F0;
    double B = (E0-Ephi+C*sin(2*phi1))/(1-cos(2*phi1));
    double A = Ephi-B*cos(2*phi1)-C*sin(2*phi1);

    double phi_min = 0.5*atan(C/B);

    // Check if we may have found a maximum rather than a minimum
    double d2E_dPhi2 = -4*B*cos(2*phi_min)-4*C*sin(2*phi_min);
    if(d2E_dPhi2 < 0) phi_min += M_PI/2;
    
    while(phi_min < -M_PI/2) phi_min += M_PI;
    while(phi_min > M_PI/2)  phi_min -= M_PI;

    if(verbose_) cout << "rotation angle = " << phi_min*360/(2*M_PI) << endl;

    //Check: Plot angular force and compare to fit
    if(doCheck)
      showGraph(phi_min,A,B,C);

    // Finally: rotate into the predicted state

    DFT_Vec y; y.zeros(N_.size());

    y.Increment_And_Scale(N_,cos(phi_min));
    y.Increment_And_Scale(T_,sin(phi_min));
    y.normalise();

    R1_.set(R0_);  R1_.Increment_And_Scale(y,dR_);
    R2_.set(R0_);  R2_.Increment_And_Scale(y,-dR_);

    DFT_Vec t(T_);
    Rotate(t,phi_min);
    T_.set(t);
    T_.normalise();

    N_.set(R1_);
    N_.Increment_And_Scale(R2_,-1);
    N_.multBy(1.0/(2*dR_));
    
    // Final check:
    if(doCheck)
      if(verbose_)
	{
	  double E;
	  cout << "Resulting angular force: " << AngularForce(0,E) << endl;
	}
  }



  
  double Energy(DFT_Vec &v) const
  {
    DFT_Vec dF; dF.zeros(v.size());
    double E;
    EnergyAndForce(v,E,dF);
    return E;
  }

  void EnergyAndForce(DFT_Vec &v, double &E, DFT_Vec &dF) const
  {
    double mu = 0.0;
    bool onlyFex = false;

    double N0 = v.size()*SMALL_VALUE*dV_;
    double sum = v.dotWith(v); 

    double a = (Nfixed_-N0)/(sum*dV_);

    DFT_Vec y(v);
    y.multBy(sqrt(a));

    density_.set_density_from_amplitude(y);

    E = dft_.calculateFreeEnergyAndDerivatives(density_, mu, dF, onlyFex);

    if(dF.size() > 0)
      {   
	double mu_eff = 0.0;
	for(long i=0;i<dF.size();i++)
	  mu_eff += dF.get(i)*v.get(i)*v.get(i);
	mu_eff *= 1.0/sum;
    
	for(int i=0;i<dF.size();i++)
	  {
	    double force = dF.get(i)-mu_eff;
	    dF.set(i,-2*a*force*v.get(i));  // minus to give the force
	  }
      }
  }

  void Rotate(DFT_Vec &v,  double theta) const
  {
    double vN = v.dotWith(N_);
    double vT = v.dotWith(T_);

    double vN_rot = vN*cos(theta)-vT*sin(theta);
    double vT_rot = vN*sin(theta)+vT*cos(theta);

    v.Increment_And_Scale(N_,vN_rot-vN);
    v.Increment_And_Scale(T_,vT_rot-vT);
  }

  double AngularForce(double dtheta, double &E) const 
  {
    DFT_Vec R1(R0_);
    R1.Increment_And_Scale(N_,dR_*cos(dtheta));
    R1.Increment_And_Scale(T_,dR_*sin(dtheta));

    DFT_Vec R2(R0_);
    R2.multBy(2);
    R2.DecrementBy(R1);
    
    double E1rot,E2rot;
    DFT_Vec dE1rot(R1_.size());
    DFT_Vec dE2rot(R2_.size());
	
    EnergyAndForce(R1,E1rot,dE1rot);
    EnergyAndForce(R2,E2rot,dE2rot);

    E = E1rot+E2rot;
    cout << "dtheta = " << dtheta << " E = " << E << endl;
    DFT_Vec ThetaStar(T_);
    Rotate(ThetaStar,dtheta);
    ThetaStar.normalise();

    DFT_Vec fPerp_Rotated(dE1rot);
    fPerp_Rotated.Increment_And_Scale(dE2rot,-1);

    // minus to give the forces (not the gradients)
    double fP = -fPerp_Rotated.dotWith(ThetaStar)*dR_;

    return fP;
  }


  bool translateDimer(double dt)
  {
    EnergyAndForce(R1_,E1_,F1_);
    EnergyAndForce(R2_,E2_,F2_);

    cout << setprecision(12) << "Before translation : E1 = " << E1_ << " E2 = " << E2_ << " av = " << (E1_+E2_)/2 << endl;
    cout << setprecision(4);
	
    double C = (F2_.dotWith(N_) - F1_.dotWith(N_))/(2*dR_);

    DFT_Vec F(F1_);
    F.IncrementBy(F2_);
    
    double FN = F.dotWith(N_);

    if(C > 0)
      {
	F.set(N_);
	F.multBy(-FN);
      } else F.Increment_And_Scale(N_,-2*FN);

    DFT_Vec dV(F);
    dV.multBy(dt);

    double dprod = V_.dotWith(dV);
    double dV2  = dV.dotWith(dV);

    V_.set(dV);
    V_.multBy(1+dprod/dV2);
    if(V_.dotWith(F) <0)V_.set(dV);    

    
    R0_.Increment_And_Scale(V_,dt);
    R1_.Increment_And_Scale(V_,dt);
    R2_.Increment_And_Scale(V_,dt);

    double E1 = Energy(R1_);
    double E2 = Energy(R2_);

    if(verbose_) 
      {
	cout << "R1_.first = " << R1_.get(0) << " R1_.second = " << R1_.get(1) << endl;
	cout << "R2_.first = " << R2_.get(0) << " R2_.second = " << R2_.get(1) << endl;
	cout << "N_[0]  = " << N_.get(0)  << " N_[1]  = " << N_.get(1) << endl;
	cout << setprecision(12);
	cout << "E1 = " << E1 << " E2 = " << E2 << " av = " << (E1+E2)/2 << " F1+F2 = " << F1_.get(0)+F2_.get(0) << "," << F1_.get(1)+F2_.get(1) << endl;
	cout << setprecision(4);
	cout << "C = " << C << " dt = " << dt <<  endl;
	density_.set(R1_);
	cout << "nAtoms = " << density_.getNumberAtoms() << endl;
      }
    double r = 0;
    for(int i=0;i<F.size();i++)
      r += fabs(F.get(i));
    
    return (r < 1e-8);

  }
  

  void showGraph(double tt, double A, double B, double C) const
  {
    {
      Grace g;
      double dd = 2*M_PI/100;
      for(int i=0;i<100;i++)
	{
	  double theta = i*dd;
	  double E;
	  double f = AngularForce(theta,E);
	  g.addPoint(theta*360/(2*M_PI),f,0); //*dR_*dR_);
	  g.addPoint(theta*360/(2*M_PI), -2*B*sin(2*theta)+2*C*cos(2*theta),1);
	  cout << "theta = " << theta << " A = " << A << " B = " << B << " C = " << C << " A+B*cos(2*theta)+C*sin(2*theta) = " << A+B*cos(2*theta)+C*sin(2*theta) << endl;
	  //	  g.addPoint(theta*360/(2*M_PI), B*cos(2*theta)+C*sin(2*theta),10);
	  g.redraw();
	}
      g.circle(0);
      g.size(0,0.2);
      g.addPoint(0,0,2);
      g.addPoint(360,0,2);
	
      double tdisplay = tt;
      if(tdisplay < 0) tdisplay = 2*M_PI+tdisplay;
      if(tdisplay > 2*M_PI) tdisplay = tdisplay - 2*M_PI;

      g.addPoint(tdisplay*360/(2*M_PI), (-2*B*sin(2*tdisplay) + 2*C*cos(2*tdisplay))/(dR_*dR_),3);
      g.square(3);
      double E;
      if(verbose_) cout << "Force = " << (-2*B*sin(2*tdisplay) + 2*C*cos(2*tdisplay))/(dR_*dR_) << " exact: " << AngularForce(tdisplay,E) << endl;
      cout << "dR_ = " << dR_ << endl;
      g.redraw();
      g.pause();
      g.close();
    }
  }

 private:
  Density &density_;
  DFT &dft_;
  
  double E0_;
  DFT_Vec R0_;
  DFT_Vec F0_;
  
  double  E1_;
  DFT_Vec R1_;
  DFT_Vec F1_;

  double  E2_;
  DFT_Vec R2_;
  DFT_Vec F2_;

  DFT_Vec V_;
  
  DFT_Vec &N_;
  DFT_Vec T_;

  double dR_;

  bool verbose_;

  double dV_;
  double Nfixed_;
};



#endif // __LUTSKO__DIMERMETHOD__
