#ifndef __LUTSKO_EXTERNAL_FIELD_
#define __LUTSKO_EXTERNAL_FIELD_

/*
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>

#include "Lattice.h"

// There is some sort of a conflict between Boost headers and <complex> so this needs to come last in order to compile.

#include <complex>


class External_Field
{
public:
  External_Field() {}
  ~External_Field() {}

  int get_species() const {return species_;}
  const DFT_Vec& get_field() const { return field_;}

  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, unsigned int version);

protected:
  DFT_Vec field_;
  int     species_ = 1;
};
*/
#endif
