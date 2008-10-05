#ifndef ARIADNE_MACROS_H
#define ARIADNE_MACROS_H

#include <sstream>
#include <stdexcept>

#define ARIADNE_THROW(except,func,msg)          \
{ \
  std::stringstream ss; \
  ss << #except " in " << func << " " << msg;    \
  throw except(ss.str()); \
} \

#define ARIADNE_ASSERT(expression) \
{ \
  bool result = (expression); \
  if(!result) { \
    ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<__PRETTY_FUNCTION__,"Assertion `" << #expression << "' failed.\n"); \
  } \
} \

#define ARIADNE_LOG(level,msg)                  \
  if(verbosity >= level) { std::clog << msg << std::flush; }

#endif // ARIADNE_MACROS_H
