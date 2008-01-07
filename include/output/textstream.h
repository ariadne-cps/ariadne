/****************************************************************************
 *            textstream.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins, Davide Bresolin
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl, bresolin@sci.univr.it
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

#ifndef ARIADNE_TEXTSTREAM_H
#define ARIADNE_TEXTSTREAM_H

/*! \file textstream.h
 *  \brief Text output.
 */
 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>

#include "base/stlio.h"
#include "numeric/conversion.h"
#include "linear_algebra/matrix.h"
#include "geometry/exceptions.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/list_set.h"


namespace Ariadne {
  namespace Output {
				
    class textstream;
    class textfstream;
    
    /*!\brief A stream for readable textual output. */
    class textstream
    {
     public:
      ~textstream();
      textstream();
      textstream(std::ostream& os);
      void redirect(std::ostream& os);
      std::ostream& ostream();
      void writenl();
 			
      template <class R> void write(const Geometry::Point<R>& pt);

     private:
      friend textstream& operator<<(textstream&, std::ostream&(*)(std::ostream&) );
     private:
      textstream(const textstream&); // no copy constructor
     private:
      std::ostream* _os_ptr;
    };


    /*!\brief A stream for readable textual output. */
    class textfstream
      : public textstream
    {
     public:
      ~textfstream();
      textfstream();
      void open(const char* fn);
      void close();
     private:
      textfstream(const textfstream&); // no copy constructor
     private:
      std::ofstream* _ofs_ptr;
    };



    inline
    std::ostream& 
    textstream::ostream()
    {
      return *this->_os_ptr; 
    }
						 
    template<class R> inline
    void 
    textstream::write(const Geometry::Point<R>& pt) 
    {
      if(pt.dimension() > 0) {
        *this->_os_ptr << approximation(pt);
      }
    }
						 

    
    inline 
    textstream& operator<<(textstream& txt, std::ostream&(*f)(std::ostream&) ) {
      std::ostream& os(*(txt._os_ptr)); os << f; return txt; 
    }
      
    template<class R> textstream& operator<<(textstream&, const Geometry::Point<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::Box<R>&);
    template<class R> textstream& operator<<(textstream&, const Geometry::GridCell<R>&);

    template<class R> textstream& operator<<(textstream&, const Geometry::PointList<R>&);
    template<class R> textstream& operator<<(textstream&, const Geometry::BoxListSet<R>&);
    template<class R> textstream& operator<<(textstream&, const Geometry::GridCellListSet<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::GridMaskSet<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::PartitionTreeSet<R>&); 

    template<class R> textstream& operator<<(textstream&, const Geometry::Zonotope<R>&);
    template<class R> textstream& operator<<(textstream&, const Geometry::Polytope<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::Polyhedron<R>&); 
    template<class BS> textstream& operator<<(textstream&, const Geometry::ListSet<BS>&); 

    template<class R> textstream& operator<<(textstream&, const Geometry::SetInterface<R>&); 

    template<class R> textstream& operator<<(textstream&, const Geometry::FiniteGrid<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::PartitionTree<R>&); 


    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Point<R>& pt) 
    { 		
      if(pt.dimension() > 0) {
        txt.write(pt);			
        txt.writenl();
      }
      return txt;
    }


    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Box<R>& bx) 
    {
      txt.ostream() << bx;
      return txt;
    }
    

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::PointList<R>& pl)
    {
      if(pl.size() > 0) {
        for (size_type i=0; i<pl.size(); i++) {
          txt.write(pl[i]);			
        }
      }
      return txt;
    }

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::BoxListSet<R>& bxls) 
    {
      txt.ostream() << bxls;
      return txt;
    }
    


 
    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::GridCell<R>& bs)
    {
      return txt << Geometry::Box<R>(bs);
    }
    

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::GridBlock<R>& bs)
    {
      return txt << Geometry::Box<R>(bs);
    }
    

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::GridCellListSet<R>& ds)
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator const_iterator;
      for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
        txt << *set_iter;
      }
      return txt;
    }
    

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::GridMaskSet<R>& ds)
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator const_iterator;
      for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
        txt << *set_iter;
      }
      return txt;
    }
    
 

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::PartitionTreeSet<R>& ds)
    {
      typedef typename Geometry::PartitionTreeSet<R>::const_iterator const_iterator;
      for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
        txt << Geometry::Box<R>(*set_iter);
      }
      return txt;
    }


    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Rectangle<R>& r)
    {
      txt.ostream() << r;
      return txt;
    }

    
    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Zonotope<R>& z)
    { 
      txt.ostream()<< z;
      return txt;
    }

       
    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Polytope<R>& p)
    {
      return txt << p.vertices();
    }

    
    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Polyhedron<R>& p)
    {
      return txt << Geometry::Polytope<Numeric::Rational>(p);
    }


    
    template<class BS> inline
    textstream&
    operator<<(textstream& txt, const Geometry::ListSet<BS>& ds)
    {
      typedef typename Geometry::ListSet<BS>::const_iterator const_iterator;
      for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
        txt << *set_iter;
      }
      return txt;
    }
    

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::SetInterface<R>& set)
    {
      return txt.ostream() << set;
    }



  }
}

#endif /* ARIADNE_TEXTSTREAM_H */
