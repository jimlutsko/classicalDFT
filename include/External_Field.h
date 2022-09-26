#ifndef __LUTSKO_EXTERNAL_FIELD_
#define __LUTSKO_EXTERNAL_FIELD_

#include "DFT_LinAlg.h"


class External_Field : public DFT_Vec
{
public:
  External_Field(): DFT_Vec() {}
  ~External_Field() {}

  int get_species() const {return 1;}
};

#endif
