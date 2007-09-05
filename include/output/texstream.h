/****************************************************************************
 *            texstream.h
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

#ifndef ARIADNE_TEXSTREAM_H
#define ARIADNE_TEXSTREAM_H

/*! \file texstream.h
 *  \brief Formatted output to TeX or LaTeX.
 */
 
#include <iosfwd>
#include <fstream>

#include "../numeric/declarations.h"
#include "../linear_algebra/declarations.h"
#include "../function/declarations.h"
#include "../geometry/declarations.h"

namespace Ariadne {
  namespace Output {
    
    class texstream;

    texstream& operator<<(texstream& txs, const char* c);
    texstream& operator<<(texstream& txs, const int& n);
    texstream& operator<<(texstream& txs, const double& x);

    template<class R> texstream& operator<<(texstream& txs, const Numeric::Integer& z);
    template<class R> texstream& operator<<(texstream& txs, const Numeric::Rational& q);
    template<class R> texstream& operator<<(texstream& txs, const Numeric::Float64& x);
    template<class R> texstream& operator<<(texstream& txs, const Numeric::FloatMP& x);

    template<class R> texstream& operator<<(texstream& txs, const Numeric::Interval<R>& ivl);
    template<class R> texstream& operator<<(texstream& txs, const LinearAlgebra::Vector<R>& v);
    template<class R> texstream& operator<<(texstream& txs, const LinearAlgebra::Matrix<R>& v);
    template<class R> texstream& operator<<(texstream& txs, const Geometry::Rectangle<R>& r);
    template<class R> texstream& operator<<(texstream& txs, const Function::Polynomial<R>& p);

    /*!\brief A stream for output to the TeX or LaTeX typesetting language. */
    class texstream {
     public:
      ~texstream();
      texstream();
      texstream(std::ostream& os);
      void redirect(std::ostream& os);
     private:
      friend texstream& operator<<(texstream&, const char&);
      friend texstream& operator<<(texstream&, const char*);
      friend texstream& operator<<(texstream&, const int&);
      friend texstream& operator<<(texstream&, const double&);
      friend texstream& operator<<(texstream&, const Numeric::Integer&);
      friend texstream& operator<<(texstream&, const Numeric::Rational&);
      friend texstream& operator<<(texstream&, const Numeric::Float64&);
      friend texstream& operator<<(texstream&, const Numeric::FloatMP&);
     private:
      std::ostream* _os_ptr;
    };

    /*!\brief A stream for file output to the TeX or LaTeX typesetting language. */
    class texfstream 
      : public texstream
    {
     public:
      texfstream();
      void open(const char* fn, const char* preamble="");
      void close();
      ~texfstream();
     private:
      std::ofstream* _ofs_ptr;
    };

  }
}

#include "texstream.template.h"

#endif /* ARIADNE_TEXSTREAM_H */


