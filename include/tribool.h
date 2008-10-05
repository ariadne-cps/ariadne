#ifndef ARIADNE_TRIBOOL_H
#define ARIADNE_TRIBOOL_H

#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>

using boost::logic::tribool;
using boost::logic::indeterminate;

inline bool definitely(tribool tb) { return tb; }
inline bool possibly(tribool tb) { return tb || indeterminate(tb); }

inline tribool operator^(tribool tb1, tribool tb2) { return (tb1&&!tb2)||(!tb1&&tb2); }

#endif /* ARIADNE_TRIBOOL_H */
