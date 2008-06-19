/****************************************************************************
 *            latexstream.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef ARIADNE_LATEXSTREAM_H
#define ARIADNE_LATEXSTREAM_H

/*! \file latexstream.h
 *  \brief Formatted output to TeX or LaTeX.
 */
 
#include <iosfwd>
#include <fstream>

#include "base/types.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"
#include "geometry/declarations.h"

namespace Ariadne {
  
    
    class latexstream;

    latexstream& operator<<(latexstream& txs, const char& c);
    latexstream& operator<<(latexstream& txs, const char* s);
    latexstream& operator<<(latexstream& txs, const int& n);
    latexstream& operator<<(latexstream& txs, const uint& n);
    latexstream& operator<<(latexstream& txs, const double& x);

    template<class R> latexstream& operator<<(latexstream& txs, const Integer& z);
    template<class R> latexstream& operator<<(latexstream& txs, const Rational& q);
    template<class R> latexstream& operator<<(latexstream& txs, const Float64& x);
    template<class R> latexstream& operator<<(latexstream& txs, const FloatMP& x);

    template<class R> latexstream& operator<<(latexstream& txs, const Interval<R>& ivl);
    template<class R> latexstream& operator<<(latexstream& txs, const Vector<R>& v);
    template<class R> latexstream& operator<<(latexstream& txs, const Matrix<R>& v);
    template<class R> latexstream& operator<<(latexstream& txs, const Box<R>& r);
    template<class R> latexstream& operator<<(latexstream& txs, const PolynomialFunction<R>& p);
    template<class R> latexstream& operator<<(latexstream& txs, const TaylorModel<R>& t);

    /*!\brief A stream for output to the TeX or LaTeX typesetting language. */
    class latexstream {
     public:
      ~latexstream();
      latexstream();
      latexstream(std::ostream& os);
      void redirect(std::ostream& os);
     private:
      friend latexstream& operator<<(latexstream&, const char&);
      friend latexstream& operator<<(latexstream&, const char*);
      friend latexstream& operator<<(latexstream&, const int&);
      friend latexstream& operator<<(latexstream&, const long int&);
      friend latexstream& operator<<(latexstream&, const unsigned int&);
      friend latexstream& operator<<(latexstream&, const unsigned long int&);
      friend latexstream& operator<<(latexstream&, const double&);
      friend latexstream& operator<<(latexstream&, const Integer&);
      friend latexstream& operator<<(latexstream&, const Rational&);
      friend latexstream& operator<<(latexstream&, const Float64&);
      friend latexstream& operator<<(latexstream&, const FloatMP&);
     private:
      std::ostream* _os_ptr;
    };

    /*!\brief A stream for file output to the TeX or LaTeX typesetting language. */
    class latexfstream 
      : public latexstream
    {
     public:
      latexfstream();
      void open(const char* fn, const char* preamble="");
      void close();
      ~latexfstream();
     private:
      std::ofstream* _ofs_ptr;
    };


} // namespace Ariadne


#include "latexstream.template.h"

#endif /* ARIADNE_LATEXSTREAM_H */


