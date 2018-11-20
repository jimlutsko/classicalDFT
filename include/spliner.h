#ifndef __SPLINER__
#define __SPLINER__

#include <vector>
#include <utility>



/**
  *  @brief UTILITY: implements the NumericalRecipes cubic spline interpolation algorithms.
  */  
class Spliner
{
 public:  
  Spliner(double const * const _x, double const * const _y, int _N1, double y0 = 1e30, double yn = 1e30);
  Spliner(Spliner &s);
Spliner() : x_(NULL), y_(NULL), y2_(NULL){}
  
  ~Spliner();

  /// calculate f
  double f(double z) const;
  /// the first derivative
  double dfdx(double z) const;
  /// the second derivative
  double d2fdx2(double z) const;
  /// the third derivative
  double d3fdx3(double z) const;

  /// integral of f
  double integrate(double a, double b);

 protected:

  void locate(double z, double &a, double &b, double &h, int &klo, int &khi) const;

  double *x_;
  double *y_;
  int N_;
  double *y2_;

};

/**
  *  @brief UTILITY: implements the "nonparametric smoothing spline" 
  */  
class smoothingSpline
{
 public:
  // x & y are the data arrays
  smoothingSpline(double const * const x, double const * const y, int N);
  ~smoothingSpline();

  // smooth using the given value of alpha
  double smooth(double alpha);

  // cross validate and return best alpha
  double crossValidate(double &alpha);

  // evaluate smoothed function
  double f(double x);

  // evaluate derivative of smoothed function
  double dfdx(double x);

 private:
  double const * const x_;
  double const * const y_;
  int N_;
  
  double *y2_;
  double *g_;


};
  

/**
  *  @brief UTILITY: 2D cubic spline interpolation implementing the NumericalRecipes 2D cubic spline interpolation algorithms. 
  */  
class Spliner2
{
 public:
  /// Grid point (i,j) has x coordinates x1[i], x2[i][j] 
  /// (so spacing of x2 can vary). The y value is y[i][j].
  Spliner2(double const * const x1, double const * const * const x2, 
	   double const * const * const y, 
	   int N1, int N2, double bc = 1e30);
  Spliner2(const std::string &file, double bc = 1e30);
  Spliner2(){}

  ~Spliner2();

  /// calculate f
  double f(double z1, double z2) const;
  double dfdx1(double z1, double z2) const;
  double dfdx2(double z1, double z2) const;
  double d2fdx1dx1(double z1, double z2) const;
  double d2fdx2dx2(double z1, double z2) const;
  double d2fdx1dx2(double z1, double z2) const;

  /// dfdx[i]
  double d1f(double z1, double z2, int i) const
    {
      if(i == 0) return dfdx1(z1,z2);
      return dfdx2(z1,z2);
    }
  /// d2fdx[i]dx[j]
  double d2f(double z1, double z2, int i, int j) const
    {
      if(i == 0 && j == 0) return d2fdx1dx1(z1,z2);
      if(i == 1 && j == 1) return d2fdx2dx2(z1,z2);
      return d2fdx1dx2(z1,z2);
    }

  double getMaxDensity();

 protected:

  void construct(double const * const x1, 
		    double const * const * const x2, 
		    double const * const * const y);


  double *x1_;
  int N1_, N2_;
  double bc_;
  double *ytmp_, *yytmp_;
  Spliner **s1_;
};

/**
  *  @brief UTILITY: implements the NumericalRecipes cubic spline interpolation algorithms. Same as Spliner but using vectors instead of real arrays
  */  
class SplinerVec
{
 public:  
  SplinerVec(const std::vector<double> & x, const std::vector<double> &y, double y0 = 1e30, double yn = 1e30);
  SplinerVec(){}
  SplinerVec(SplinerVec const &s);
  ~SplinerVec();

  void initialize(const std::vector<double> & x, const std::vector<double> &y, double y0 = 1e30, double yn = 1e30);

  /// calculate f
  double f(double z) const;
  /// the first derivative
  double dfdx(double z) const;
  /// the second derivative
  double d2fdx2(double z) const;
  /// the third derivative
  double d3fdx3(double z) const;

  /// integral of f
  double integrate(double a, double b);

  double getXMax() const { return x_.back();}
  double getXMin() const { return x_.front();}

 protected:

  void locate(double z, double &a, double &b, double &h, int &klo, int &khi) const;

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> y2_;
};

/**
  *  @brief UTILITY: implements the NumericalRecipes cubic spline interpolation algorithms. Same as Spliner but using the pair template object.
  */  
class SplinerPair
{
public:  
    SplinerPair(const std::vector< std::pair<double,double> > &data, double y0 = 1e30, double yn = 1e30)
    {initialize(data,y0,yn);}
    SplinerPair(){}
    SplinerPair(SplinerPair const &s);
    ~SplinerPair();

    void initialize(std::vector< std::pair<double,double> > const &data, double y0 = 1e30, double yn = 1e30);

    /// calculate f
    double f(double z) const;
    /// the first derivative
    double dfdx(double z) const;
    /// the second derivative
    double d2fdx2(double z) const;
    /// the third derivative
    double d3fdx3(double z) const;

    /// integral of f
    double integrate(double a, double b);

protected:

    double x(int i) const { return data_[i].first;}
    double y(int i) const { return data_[i].second;}

    void locate(double z, double &a, double &b, double &h, int &klo, int &khi) const;

    std::vector< std::pair<double,double> > data_;
    std::vector<double> y2_;
};

/**
  *  @brief UTILITY: implements the NumericalRecipes cubic spline interpolation algorithms for 2D splines. Same as Spliner2 but it allows for different numbers of points for each value of x.
  */  
class Spliner2General
{
public:
    /// Grid point (i,j) has x coordinates x1[i], x2[i][j] 
    /// (so spacing of x2 can vary). The y value is y[i][j].
    Spliner2General( std::vector< std::pair < double, std::vector< std::pair<double,double> > > > const &data, double bc = 1e30,
		     double bc20 = 1e30, double bc21 = 1e30);
    Spliner2General(){}
    void initialize(std::vector< std::pair< double, std::vector< std::pair<double,double> > > >  const &data,double bc = 1e30, double bc1 = 1e30, double bc2 = 1e30);

    ~Spliner2General(){}

    /// calculate f
    double f(double z1, double z2) const;
    double dfdx1(double z1, double z2) const;
    double dfdx2(double z1, double z2) const;
    double d2fdx1dx1(double z1, double z2) const;
    double d2fdx2dx2(double z1, double z2) const;
    double d2fdx1dx2(double z1, double z2) const;

    /// dfdx[i]
    double d1f(double z1, double z2, int i) const
    {
	if(i == 0) return dfdx1(z1,z2);
	return dfdx2(z1,z2);
    }
    /// d2fdx[i]dx[j]
    double d2f(double z1, double z2, int i, int j) const
    {
	if(i == 0 && j == 0) return d2fdx1dx1(z1,z2);
	if(i == 1 && j == 1) return d2fdx2dx2(z1,z2);
	return d2fdx1dx2(z1,z2);
    }

    double getMaxDensity();

protected:
    std::vector< std::pair< double, std::vector< std::pair<double,double> > > > data_;

    double bc_;
    std::vector< SplinerPair > s1_;
};

#endif //__SPLINER__

