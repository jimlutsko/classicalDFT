#include "dft_lib/exceptions/parameter_exceptions.h"

#include <utility>

namespace dft_core
{
  namespace exception
  {
    NegativeParameterException::NegativeParameterException(std::string msg)
    : WrongParameterException(std::move(msg))
    { }
  }
}