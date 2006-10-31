/***************************************************************************
 *            tribool.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
 *
 ****************************************************************************/

#ifndef _ARIADNE_TRIBOOL_H
#define _ARIADNE_TRIBOOL_H

#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>

namespace Ariadne {
  namespace Base {
    using boost::logic::tribool;
    using boost::logic::indeterminate;

    inline bool possibly(tribool tb) { return tb || indeterminate(tb); }
  }
}

#endif
