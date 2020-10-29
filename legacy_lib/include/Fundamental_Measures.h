#ifndef __LUTSKO_FUNDAMENTAL_MEASURES__
#define __LUTSKO_FUNDAMENTAL_MEASURES__


class FundamentalMeasures
{
 public:
  FundamentalMeasures()
    {
      eta = s0 = s1 = s2 = 0.0;
      v1[0] = v1[1] = v1[2] = 0.0;
      v2[0] = v2[1] = v2[2] = 0.0;
      T[0][0] = T[0][1] = T[0][2] = 0.0;
      T[1][0] = T[1][1] = T[1][2] = 0.0;
      T[2][0] = T[2][1] = T[2][2] = 0.0;
      calculate_derived_quantities();      
    }
  FundamentalMeasures(double density, double hsd) {fillUniform(density,hsd);}

  void fillUniform(double density, double hsd = 1)
  {
      eta = M_PI*density*hsd*hsd*hsd/6;
      s0  = M_PI*density;
      s1 = s0*hsd;
      s2 = s1*hsd;
      v1[0] = v1[1] = v1[2] = 0.0;
      v2[0] = v2[1] = v2[2] = 0.0;
      T[0][0] = T[1][1] = T[2][2] = s2/3;
      T[0][1] = T[0][2] = T[1][0] = T[1][2] = T[2][0] = T[2][1] = 0.0;
      calculate_derived_quantities();
    }


  FundamentalMeasures(const FundamentalMeasures& b)
    {
      eta = b.eta;
      s0 = b.s0;
      s1 = b.s1;
      s2 = b.s2;
      for(int i=0;i<3;i++)
	{
	  v1[i] = b.v1[i];
	  v2[i] = b.v2[i];
	  vT[i] = b.vT[i];       
	  for(int j=0;j<3;j++)
	    {
	      T[i][j] = b.T[i][j];
	      TT[i][j] = b.TT[i][j];
	    }
	}
      v1_v2 = b.v1_v2;
      v2_v2 = b.v2_v2;
      vTv = b.vTv;
      T2 = b.T2;
      T3 = b.T3;
    }
  
  FundamentalMeasures(double f[]){fill(f);}  

  void calculate_derived_quantities()
  {
    vT[0] = T[0][0]*v2[0] +T[0][1]*v2[1] +T[0][2]*v2[2];
    vT[1] = T[1][0]*v2[0] +T[1][1]*v2[1] +T[1][2]*v2[2];
    vT[2] = T[2][0]*v2[0] +T[2][1]*v2[1] +T[2][2]*v2[2];

    TT[0][0] = T[0][0]*T[0][0] +T[0][1]*T[1][0] +T[0][2]*T[2][0];
    TT[1][0] = T[1][0]*T[0][0] +T[1][1]*T[1][0] +T[1][2]*T[2][0];
    TT[2][0] = T[2][0]*T[0][0] +T[2][1]*T[1][0] +T[2][2]*T[2][0];

    TT[0][1] = T[0][0]*T[0][1] +T[0][1]*T[1][1] +T[0][2]*T[2][1];
    TT[1][1] = T[1][0]*T[0][1] +T[1][1]*T[1][1] +T[1][2]*T[2][1];
    TT[2][1] = T[2][0]*T[0][1] +T[2][1]*T[1][1] +T[2][2]*T[2][1];

    TT[0][2] = T[0][0]*T[0][2] +T[0][1]*T[1][2] +T[0][2]*T[2][2];
    TT[1][2] = T[1][0]*T[0][2] +T[1][1]*T[1][2] +T[1][2]*T[2][2];
    TT[2][2] = T[2][0]*T[0][2] +T[2][1]*T[1][2] +T[2][2]*T[2][2];

    v1_v2  = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]; 
    v2_v2  = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]; 
    vTv    = v2[0]*vT[0]+v2[1]*vT[1]+v2[2]*vT[2]; 
    T2     = TT[0][0]+TT[1][1]+TT[2][2]; 
    T3     = TT[0][0]*T[0][0]+TT[0][1]*T[1][0]+TT[0][2]*T[2][0]
      +TT[1][0]*T[0][1]+TT[1][1]*T[1][1]+TT[1][2]*T[2][1]
      +TT[2][0]*T[0][2]+TT[2][1]*T[1][2]+TT[2][2]*T[2][2];
    
  }

  void fill(double f[])
  {
    eta   = f[0];
    s0    = f[1];
    s1    = f[2];
    s2    = f[3];
    v1[0] = f[4];
    v1[1] = f[5];
    v1[2] = f[6];
    v2[0] = f[7];
    v2[1] = f[8];
    v2[2] = f[9];

    T[0][0] = f[10];
    T[0][1] = f[11];
    T[0][2] = f[12];

    T[1][0] = f[13];
    T[1][1] = f[14];
    T[1][2] = f[15];

    T[2][0] = f[16];
    T[2][1] = f[17];
    T[2][2] = f[18];  

    calculate_derived_quantities();  
  }
  
  //  void copy_to(double f[])
  vector<double> get_as_vector()
  {
    vector<double> f;
    f.push_back(eta);
    f.push_back(s0);
    f.push_back(s1);
    f.push_back(s2);
    f.push_back(v1[0]);
    f.push_back(v1[1]);
    f.push_back(v1[2]);
    f.push_back(v2[0]);
    f.push_back(v2[1]);
    f.push_back(v2[2]);

    f.push_back(T[0][0]);
    f.push_back(T[0][1]);
    f.push_back(T[0][2]);
    
    f.push_back(T[1][0]);
    f.push_back(T[1][1]);
    f.push_back(T[1][2]);

    f.push_back(T[2][0]);
    f.push_back(T[2][1]);
    f.push_back(T[2][2]);

    return f;
  }

  //  void fillHomogeneousWeights(double d)
  //  {
  //  }


  // fundamental measures
  double eta = 0.0;
  double s0 = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  double v1[3] = {0.0,0.0,0.0};
  double v2[3] = {0.0,0.0,0.0};
  double T[3][3] = {{0.0,0.0,0.0},
		    {0.0,0.0,0.0},
		    {0.0,0.0,0.0}};

  // some useful derived quantities
  double v1_v2 = 0.0; ///< v1 dot v2
  double v2_v2 = 0.0; ///< v2 dot v2
  double vTv   = 0.0; ///< v dot T dot v
  double T2    = 0.0; ///< Tr(T^2)
  double T3    = 0.0; ///< Tr(T^3)

  double vT[3] = {0.0,0.0,0.0};

  
  double TT[3][3] = {{0.0,0.0,0.0},
		     {0.0,0.0,0.0},
		     {0.0,0.0,0.0}};  
};
  










#endif
