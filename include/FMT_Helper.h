class FMT_Helper
{
 public:
  FMT_Helper(){}
  ~FMT_Helper(){}

  static double G_eta(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
  {
    double Iy[5];
    double Iz[5];
    double Jy[5];
    double Jz[5];

    // A temporary workaround to avoid NAN.
    R += 1e-15;
  
    getI(X, sqrt(R*R-Vy*Vy), Iy);
    getI(X, sqrt(R*R-Vz*Vz), Iz);
    getJ(X,Vy,R,Jy);
    getJ(X,Vz,R,Jz);

    double A = Ty*Tz*Vy*Vz+0.25*Vy*Vy*Vz*Vz+0.125*R*R*R*R
      -(1.0/6)*(Ty*Vy*Vy*Vy+Tz*Vz*Vz*Vz) -0.5*Tz*Vz*Vy*Vy-0.5*Ty*Vy*Vz*Vz
      +0.125*Vy*Vy*Vy*Vy+0.125*Vz*Vz*Vz*Vz+0.5*R*R*(Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*M_PI*Ty*Tz);
    double B = 0.5*Ty*Vy+0.5*Tz*Vz-0.25*Vy*Vy-0.25*Vz*Vz+(M_PI/4)*Ty*Tz+0.25*R*R;
    double C = 0.5*Ty*Vy+(1.0/3)*(R*R-Vy*Vy);
    double D = 0.5*Tz*Vz+(1.0/3)*(R*R-Vz*Vz);

    double g = Tx*A*X-0.5*A*X*X-(1.0/3)*B*Tx*X*X*X+0.25*X*X*X*X*B
      +0.025*Tx*X*X*X*X*X-(1.0/48)*X*X*X*X*X*X;

    if(isnan(g))
      throw std::runtime_error("Here 1");
  
    g -= 0.5*Ty*Tz*(Tx*R*R*Jy[0]-R*R*Jy[1]-Tx*Jy[2]+Jy[3]);
    if(isnan(g))
      throw std::runtime_error("Here 2");  
    g -= 0.5*Ty*Tz*(Tx*R*R*Jz[0]-R*R*Jz[1]-Tx*Jz[2]+Jz[3]);
    if(isnan(g))
      {
	cout << Jz[0] << " " << Jz[1] << " " << Jz[2] << " " << Jz[3] << endl;
	throw std::runtime_error("Here 3");
      }
    g -= (Tx*Tz*C*Iy[0]-Tz*C*Iy[1]-(1.0/3)*Tx*Tz*Iy[2]+(1.0/3)*Tz*Iy[3]);
    if(isnan(g))
      throw std::runtime_error("Here 4");  
    g -= (Tx*Ty*D*Iz[0]-Ty*D*Iz[1]-(1.0/3)*Tx*Ty*Iz[2]+(1.0/3)*Ty*Iz[3]);

    return g;
  }

  static double G_s(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
  {
    double Iy[5];
    double Iz[5];
    double Jy[5];
    double Jz[5];

    // A temporary workaround to avoid NAN.
    R += 1e-15;
  
    getI(X, sqrt(R*R-Vy*Vy), Iy);
    getI(X, sqrt(R*R-Vz*Vz), Iz);
    getJ(X,Vy,R,Jy);
    getJ(X,Vz,R,Jz);

    double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
    double g = 0.125*R*X*X*X*X-(1.0/6)*R*Tx*X*X*X-0.5*R*A*X*X+R*A*Tx*X;

    g -= R*Tz*(Tx*Iy[0]-Iy[1]);
    g -= R*Ty*(Tx*Iz[0]-Iz[1]);

    g -= R*Ty*Tz*(Tx*Jy[0]-Jy[1]);
    g -= R*Ty*Tz*(Tx*Jz[0]-Jz[1]);

    return g;
  }

  static double G_vx(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
  {
    double Iy[5];
    double Iz[5];
    double Jy[5];
    double Jz[5];

    // A temporary workaround to avoid NAN.
    R += 1e-15;  

    getI(X, sqrt(R*R-Vy*Vy), Iy);
    getI(X, sqrt(R*R-Vz*Vz), Iz);
    getJ(X,Vy,R,Jy);
    getJ(X,Vz,R,Jz);

    double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
    double g = 0.1*R*X*X*X*X*X-0.125*R*Tx*X*X*X*X-(1.0/3)*R*A*X*X*X+0.5*R*A*Tx*X*X;

    g -= R*Tz*(Tx*Iy[1]-Iy[2]);
    g -= R*Ty*(Tx*Iz[1]-Iz[2]);

    g -= R*Ty*Tz*(Tx*Jy[1]-Jy[2]);
    g -= R*Ty*Tz*(Tx*Jz[1]-Jz[2]);

  
    return g;
  }

  static double G_txx(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
  {
    double Iy[5];
    double Iz[5];
    double Jy[5];
    double Jz[5];

    // A temporary workaround to avoid NAN.
    R += 1e-15;
  
    getI(X, sqrt(R*R-Vy*Vy), Iy);
    getI(X, sqrt(R*R-Vz*Vz), Iz);
    getJ(X,Vy,R,Jy);
    getJ(X,Vz,R,Jz);

    double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
    double g = (1.0/12)*R*X*X*X*X*X*X-0.1*R*Tx*X*X*X*X*X-0.25*R*A*X*X*X*X+(1.0/3)*R*A*Tx*X*X*X;

    g -= R*Tz*(Tx*Iy[2]-Iy[3]);
    g -= R*Ty*(Tx*Iz[2]-Iz[3]);

    g -= R*Ty*Tz*(Tx*Jy[2]-Jy[3]);
    g -= R*Ty*Tz*(Tx*Jz[2]-Jz[3]);

    return g;
  }

  static double G_txy(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
  {
    double Iy[5];
    double Iz[5];
    double Jy[5];
    double Jz[5];

    // A temporary workaround to avoid NAN.
    R += 1e-15;
  
    getI(X, sqrt(R*R-Vy*Vy), Iy);
    getI(X, sqrt(R*R-Vz*Vz), Iz);
    getJ(X,Vy,R,Jy);
    getJ(X,Vz,R,Jz);

    double A = (1.0/6)*(3*Ty*R*R+2*Vy*Vy*Vy-3*Ty*Vy*Vy-3*Ty*Vz*Vz+6*Ty*Tz*Vz+1.5*M_PI*R*R*Tz);
    double B = 0.25*M_PI*Tz+0.5*Ty;
    double C = -(1.0/6)*(2*R*R-2*Vz*Vz+3*Tz*Vz);
  
    double g = -0.5*R*A*Tx*X*X+(1.0/3)*R*A*X*X*X+0.25*R*B*Tx*X*X*X*X-0.2*R*B*X*X*X*X*X;

    g += 0.5*R*Tz*(Tx*R*R*Jy[1]-R*R*Jy[2]-Tx*Jy[3]+Jy[4]);
    g += 0.5*R*Tz*(Tx*R*R*Jz[1]-R*R*Jz[2]-Tx*Jz[3]+Jz[4]);

    g += 0.5*R*Tz*(2*Ty-Vy)*(Tx*Iy[1]-Iy[2]);
    g += R*(-Tx*C*Iz[1]+C*Iz[2]-(1.0/3)*Tx*Iz[3]+(1.0/3)*Iz[4]);

    return g;
  }  

  static double G_I1x(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
  {
    double Iy[5];
    double Iz[5];
    double Jy[5];
    double Jz[5];

    // A temporary workaround to avoid NAN.
    R += 1e-15;
  
    getI(X, sqrt(R*R-Vy*Vy), Iy);
    getI(X, sqrt(R*R-Vz*Vz), Iz);
    getJ(X,Vy,R,Jy);
    getJ(X,Vz,R,Jz);

    double A = Ty*Tz*Vy*Vz+0.25*Vy*Vy*Vz*Vz+0.125*R*R*R*R
      -(1.0/6)*(Ty*Vy*Vy*Vy+Tz*Vz*Vz*Vz) -0.5*Tz*Vz*Vy*Vy-0.5*Ty*Vy*Vz*Vz
      +0.125*Vy*Vy*Vy*Vy+0.125*Vz*Vz*Vz*Vz+0.5*R*R*(Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*M_PI*Ty*Tz);
    double B = 0.5*Ty*Vy+0.5*Tz*Vz-0.25*Vy*Vy-0.25*Vz*Vz+(M_PI/4)*Ty*Tz+0.25*R*R;
    double C = 0.5*Ty*Vy+(1.0/3)*(R*R-Vy*Vy);
    double D = 0.5*Tz*Vz+(1.0/3)*(R*R-Vz*Vz);

    double X4 = X*X*X*X;
    
    double g = 0.5*Tx*A*X*X-(1.0/3)*A*X*X*X-0.25*B*Tx*X4+0.2*B*X4*X
      +(1.0/48)*Tx*X4*X*X-(1.0/56)*X4*X*X*X;

    if(isnan(g))
      throw std::runtime_error("Here 1");
  
    g -= 0.5*Ty*Tz*(Tx*R*R*Jy[1]-R*R*Jy[2]-Tx*Jy[3]+Jy[4]);
    if(isnan(g))
      throw std::runtime_error("Here 2");  
    g -= 0.5*Ty*Tz*(Tx*R*R*Jz[1]-R*R*Jz[2]-Tx*Jz[3]+Jz[4]);
    if(isnan(g))
      {
	cout << Jz[0] << " " << Jz[1] << " " << Jz[2] << " " << Jz[3] << endl;
	throw std::runtime_error("Here 3");
      }
    g -= (Tx*Tz*C*Iy[1]-Tz*C*Iy[2]-(1.0/3)*Tx*Tz*Iy[3]+(1.0/3)*Tz*Iy[4]);
    if(isnan(g))
      throw std::runtime_error("Here 4");  
    g -= (Tx*Ty*D*Iz[1]-Ty*D*Iz[2]-(1.0/3)*Tx*Ty*Iz[3]+(1.0/3)*Ty*Iz[4]);
  
    return g;
  }

  
  static double G_I2(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
  {
    double Iy[7];
    double Iz[7];
    double Jy[7];
    double Jz[7];

    // A temporary workaround to avoid NAN.
    R += 1e-15;

    bool bExtended = true;
    
    getI(X, sqrt(R*R-Vy*Vy), Iy, bExtended);
    getI(X, sqrt(R*R-Vz*Vz), Iz, bExtended);
    getJ(X,Vy,R,Jy, bExtended);
    getJ(X,Vz,R,Jz, bExtended);

    double A = Ty*Tz*Vy*Vz+0.25*Vy*Vy*Vz*Vz+0.125*R*R*R*R
      -(1.0/6)*(Ty*Vy*Vy*Vy+Tz*Vz*Vz*Vz) -0.5*Tz*Vz*Vy*Vy-0.5*Ty*Vy*Vz*Vz
      +0.125*Vy*Vy*Vy*Vy+0.125*Vz*Vz*Vz*Vz+0.5*R*R*(Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*M_PI*Ty*Tz);
    double B = 0.5*Ty*Vy+0.5*Tz*Vz-0.25*Vy*Vy-0.25*Vz*Vz+(M_PI/4)*Ty*Tz+0.25*R*R;
    double C = 0.5*Ty*Vy+(1.0/3)*(R*R-Vy*Vy);
    double D = 0.5*Tz*Vz+(1.0/3)*(R*R-Vz*Vz);

    double X4 = X*X*X*X;
    double g = (1.0/3)*Tx*A*X*X*X-0.25*A*X4-0.2*B*Tx*X4*X+(1.0/6)*X4*X*X*B
      +(1.0/56)*Tx*X4*X*X*X-0.015625*X4*X4;

    if(isnan(g))
      throw std::runtime_error("Here 1");
  
    g -= 0.5*Ty*Tz*(Tx*R*R*Jy[2]-R*R*Jy[3]-Tx*Jy[4]+Jy[5]);
    if(isnan(g))
      throw std::runtime_error("Here 2");  
    g -= 0.5*Ty*Tz*(Tx*R*R*Jz[2]-R*R*Jz[3]-Tx*Jz[4]+Jz[5]);
    if(isnan(g))
      {
	cout << Jz[0] << " " << Jz[1] << " " << Jz[2] << " " << Jz[3] << " " << Jz[4] << " " << Jz[5] << endl;
	throw std::runtime_error("Here 3");
      }
    g -= (Tx*Tz*C*Iy[2]-Tz*C*Iy[3]-(1.0/3)*Tx*Tz*Iy[4]+(1.0/3)*Tz*Iy[5]);
    if(isnan(g))
      throw std::runtime_error("Here 4");  
    g -= (Tx*Ty*D*Iz[2]-Ty*D*Iz[3]-(1.0/3)*Tx*Ty*Iz[4]+(1.0/3)*Ty*Iz[5]);

    return g;
  }



  
  static double SMALL_NUM;    
  
 private:

  static void getI(double X, double A, double I[], bool bExtended = false)
  {
    double a = asin(X/A);
    double b = sqrt(fabs(A*A-X*X));
  
    I[0] = 0.5*X*b+0.5*A*A*a;
    I[1] = -(1.0/3)*(A*A-X*X)*b;
    I[2] = 0.125*X*(2*X*X-A*A)*b+0.125*A*A*A*A*a;
    I[3] = -(1.0/15)*(A*A-X*X)*(3*X*X+2*A*A)*b;
    I[4] = (1.0/48)*(8*X*X*X*X-2*A*A*X*X-3*A*A*A*A)*X*b+0.0625*A*A*A*A*A*A*a;

    if(bExtended)
      {
        I[5] = -(1.0/105)*(A*A-X*X)*(15*X*X*X*X+12*A*A*X*X+8*A*A*A*A)*b;
        I[6] = (1.0/384)*(48*X*X*X*X*X*X-8*A*A*X*X*X*X-10*A*A*A*A*X*X-15*A*A*A*A*A*A)*X*b+0.0390625*A*A*A*A*A*A*A*A*a;
      }
  }

  static void getJ(double X, double V, double R, double J[], bool bExtended = false)
  {
    double aV  = asin(max(-1.0,min(1.0,V/max(SMALL_NUM,sqrt(R*R-X*X)))));
    double aX  = asin(max(-1.0,min(1.0,X/sqrt(R*R-V*V))));
    double aVX = asin(max(-1.0,min(1.0,X*V/max(SMALL_NUM,sqrt((R*R-V*V)*(R*R-X*X))))));
    double b   = sqrt(fabs(R*R-V*V-X*X));
  
    J[0] = X*aV+V*aX-R*aVX;
    J[1] = -0.5*V*b+0.5*(X*X-R*R)*aV;
    J[2] = -(1.0/6)*V*(V*V-3*R*R)*aX+(1.0/3)*X*X*X*aV-(1.0/3)*R*R*R*aVX-(1.0/6)*V*X*b;
    J[3] = 0.25*(X*X*X*X-R*R*R*R)*aV+(M_PI/8)*R*R*R*R-(V/12.0)*(5*R*R+X*X-2*V*V)*b;
    J[4] = 0.025*X*V*(3*V*V-2*X*X-7*R*R)*b+0.025*V*(3*V*V*V*V-10*V*V*R*R+15*R*R*R*R)*aX+0.2*X*X*X*X*X*aV-0.2*R*R*R*R*R*aVX;

    if(bExtended)
      {
	J[5] = (1.0/6)*(X*X*X*X*X*X-R*R*R*R*R*R)*aV-(V/90.0)*(33*R*R*R*R-26*R*R*V*V+9*R*R*X*X+8*V*V*V*V-4*V*V*X*X+3*X*X*X*X)*b;
	J[6] = -(1.0/336)*X*V*(15*V*V*V*V+8*X*X*X*X+57*R*R*R*R-10*V*V*X*X+22*R*R*X*X-48*R*R*V*V)*b+0.0625*V*(R*R-V*V)*(5*R*R*R*R+3*V*V*V*V)*aX
	  +(1.0/7)*V*V*V*V*V*V*V*aX+(1.0/7)*X*X*X*X*X*X*X*aV-(1.0/7)*R*R*R*R*R*R*R*aVX;
      }
    
    if(isnan(J[1]) || isnan(J[2]) || isnan(J[3]) || isnan(J[4]))
      {
	cout << X << " " << V << " " << R << " | " << b << " " << R*R-V*V << " " << R*R-X*X << " " << R*R - V*V -X*X << endl;
	throw std::runtime_error("done");
      }
  
  }










};
