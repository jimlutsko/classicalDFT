#ifndef CLASSICALDFT_GRACE_EXCEPTION_H
#define CLASSICALDFT_GRACE_EXCEPTION_H

#include <string>
#include <utility>

namespace dft_core
{
  namespace grace_plot
  {

    /**
     * @brief Exceptions related with the Grace wrapper of `xmgrace`
     *
     * This class will serve as the base class for the exceptions we want to register regarding undesired or uncontrolled
     * behaviours when dealing with `xmgrace`
     */
    class GraceException: public std::exception
    {
      private:
        /// The error message which will be deliver when the exception is trhown
        std::string error_message_;

      protected:
        /// A setter of the error message for internal use in constructors of derived classes
        void set_error_message(const std::string& message) { error_message_ = message; }

      public:
        /// Explicit constructor of the exception, to avoid unintentional implicit conversions
        explicit GraceException(std::string msg): error_message_(std::move(msg)) {}

        /// Standard default destructor to be used, it should not throw an exception and if it does the program just crash
        ~GraceException() noexcept override = default;

        /// Overrides the `what` method which is used when catching the message of a `std::exception`
        virtual const char* what() const noexcept override { return error_message_.c_str(); }

        /// Gets the error message as a constant reference to the private field `error_message_`
        const std::string& error_message() { return error_message_; }
    };

    /**
     * @brief Exception to be thrown when detecting that `xmgrace` communication is not opened (`GraceIsOpen() == false`)
     *
     * This class is just a derived class from the general `GraceException` which particularises the message to
     * specify that there is no connection established with `xmgrace` just yet.
     */
    class GraceNotOpenedException: public GraceException
    {
      public:
        GraceNotOpenedException();
    };

    /**
     * @brief Exception to be thrown when detecting that `xmgrace` communication fails
     *
     * This class is just a derived class from the general `GraceException` which particularises the message to
     * specify that there is a communication error when trying to open the communication pipe or when sending a
     * command to an open communication.
     */
    class GraceCommunicationFailedException: public GraceException
    {
      public:
        GraceCommunicationFailedException();
        explicit GraceCommunicationFailedException(const std::string& msg);
    };

  }
}

#endif
