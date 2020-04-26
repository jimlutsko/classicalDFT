#include "dft_lib/exceptions/grace_exception.h"

namespace dft_core
{
  namespace grace_plot
  {
    GraceNotOpenedException::GraceNotOpenedException() : GraceException(std::string())
    {
      set_error_message("No grace subprocess currently connected!");
    }

    GraceCommunicationFailedException::GraceCommunicationFailedException()
    : GraceException(std::string())
    {
      set_error_message("There was a problem while communicating with Grace.");
    }

    GraceCommunicationFailedException::GraceCommunicationFailedException(const std::string& msg)
    : GraceException(msg) {}
  }
}

