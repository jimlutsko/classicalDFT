#ifndef __LUTSKO_DISPLAY__
#define  __LUTSKO_DISPLAY__


#include <mgl2/mgl.h>

/**
  *  @brief UTILITY: this class can be specialized to display Density objects using external libraries. This is currently implemented using mgl. 
  */  
class Display
{
 public:
  Display(unsigned nx, unsigned ny){ gr_ = new mglGraph(); data_2D_.Create(nx,ny);}
  ~Display(){if(gr_) delete gr_;}

  void setData(unsigned n, double d) { data_2D_[n] = d;}

  
  void doDisplay(string &title, string &file)
  {
    if(!gr_) return;

    //clear the window
    gr_->Clf();	

    // basic window formating
    gr_->Box();  
    gr_->Alpha(false);
    gr_->SetRange('c', 0.0, 1.0);

    // set the color scheme
    //	  gr_->Dens(a,"kRryw");
    //	  gr_->Dens(a,"UbcyqR");
    gr_->Dens(data_2D_,"kBbcw");

    /*
    // Write a title
    char str[48];	
    snprintf(str,48,"calls = %d F = %lf N = %lf",calls_, F_, density_.getNumberAtoms());
    */
    gr_->Puts(mglPoint(0,1.1),title.c_str());

    
    gr_->WriteFrame(file.c_str());
  }

 protected:
  mglGraph *gr_;
  mglData data_2D_;

};


#endif //  __LUTSKO_DISPLAY__
