#include <exception>
#include <stdexcept>


class NotImplemented : public std::logic_error {
 public:
  NotImplemented(const std::string& str) : std::logic_error(str) { }
};

class IncompatibleSizes : public std::runtime_error {
 public:
  IncompatibleSizes(const std::string& str) : std::runtime_error(str) { }
};

