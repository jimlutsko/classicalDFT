#ifndef __Lutsko_TimeStamp__
#define __Lutsko_TimeStamp__


#include<stdio.h>
#include<time.h>



class TimeStamp
{
 public:
  TimeStamp(){}
  ~TimeStamp(){}

  /// Allow printing of the description
  friend std::ostream& operator<< (std::ostream& o, const TimeStamp& ps)
    {
      ps.printOn(o);
      return o;
    }

 protected:

  /// A unique description of this class
  virtual void printOn(std::ostream& o) const
  {
    char timestamp[100];
    time_t mytime;
    struct tm *mytm;
    mytime=time(NULL);
    mytm=localtime(&mytime);
    
    strftime(timestamp,sizeof timestamp,"%a, %d %b %Y %H:%M:%S %z",mytm);
    o << timestamp;
  }

};


#endif // __Lutsko_TimeStamp__
