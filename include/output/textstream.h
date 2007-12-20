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
#include "geometry/rectangle.h"
#include "geometry/rectangular_set.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedral_set.h"
#include "geometry/partition_tree_set.h"
#include "system/affine_map.h"
 


namespace Ariadne {
  namespace Output {
				
    class textstream;
    class textfstream;
    
    template<class BS> std::string summary(const Geometry::ListSet<BS>& ls);
    template<class R> std::string summary(const Geometry::GridCellListSet<R>& gcls);
    template<class R> std::string summary(const Geometry::GridMaskSet<R>& gms);



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
    template<class R> textstream& operator<<(textstream&, const Geometry::PointList<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::Rectangle<R>&);
    template<class R, class Tag> textstream& operator<<(textstream&, const Geometry::Zonotope<R,Tag>&);
    template<class R> textstream& operator<<(textstream&, const Geometry::Polytope<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::Polyhedron<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::RectangularSet<R>&);
    template<class R> textstream& operator<<(textstream&, const Geometry::PolyhedralSet<R>&);
    template<class BS> textstream& operator<<(textstream&, const Geometry::ListSet<BS>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::GridCellListSet<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::GridMaskSet<R>&); 
    template<class R> textstream& operator<<(textstream&, const Geometry::PartitionTreeSet<R>&); 
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
    operator<<(textstream& txt, const Geometry::PointList<R>& pl)
    {
      if(pl.size() > 0) {
        for (size_type i=0; i<pl.size(); i++) {
          txt.write(pl[i]);			
        }
        txt.writenl();
      }
      return txt;
    }


    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Rectangle<R>& r) 
    {
      return txt << r.vertices();
    }
    

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::RectangularSet<R>& rs)
    {
      return txt << Geometry::Rectangle<R>(rs);
    }

    
    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Zonotope< Numeric::Interval<R>,Numeric::Interval<R> >& iz)
    { 
      return txt << iz.vertices();
    }

       
    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Zonotope<Numeric::Interval<R>,R>& ez)
    { 
      return txt << ez.vertices();
    }
      
 
    template<class R, class Tag> inline
    textstream&
    operator<<(textstream& txt, const Geometry::Zonotope<R,Tag>& z)
    {  
      txt.ostream() << z;
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


    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::PolyhedralSet<R>& ps)
    {
      return txt << Geometry::Polyhedron<R>(ps);
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
    operator<<(textstream& txt, const Geometry::GridCell<R>& bs)
    {
      return txt << Geometry::Rectangle<R>(bs);
    }
    

    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::GridBlock<R>& bs)
    {
      return txt << Geometry::Rectangle<R>(bs);
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
        txt << Geometry::Rectangle<R>(*set_iter);
      }
      return txt;
    }


    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::SetInterface<R>& set)
    {
      using namespace Geometry;
      typedef Numeric::Interval<R> I;
      
      if(dynamic_cast<const RectangularSet<R>*>(&set)) {
        return txt << dynamic_cast<const RectangularSet<R>&>(set);
      } else if(dynamic_cast<const PolyhedralSet<R>*>(&set)) {
        return txt << dynamic_cast<const PolyhedralSet<R>&>(set);
      } else if(dynamic_cast<const ListSet< Rectangle<R> >*>(&set)) {
        return txt << dynamic_cast<const ListSet< Rectangle<R> >&>(set);
      } else if(dynamic_cast<const ListSet< Zonotope<R,R> >*>(&set)) {
        return txt << dynamic_cast<const ListSet< Zonotope<R,R> >&>(set);
      } else if(dynamic_cast<const ListSet< Zonotope<I,R> >*>(&set)) {
        return txt << dynamic_cast<const ListSet< Zonotope<I,R> >&>(set);
      } else if(dynamic_cast<const ListSet< Zonotope<I,I> >*>(&set)) {
        return txt << dynamic_cast<const ListSet< Zonotope<I,I> >&>(set);
      } else if(dynamic_cast<const GridCellListSet<R>*>(&set)) {
        return txt << dynamic_cast<const GridCellListSet<R>&>(set);
      } else if(dynamic_cast<const GridMaskSet<R>*>(&set)) {
        return txt << dynamic_cast<const GridMaskSet<R>&>(set);
      } else if(dynamic_cast<const PartitionTreeSet<R>*>(&set)) {
        return txt << dynamic_cast<const PartitionTreeSet<R>&>(set);
      }  else {
        Rectangle<R> bb;
        try {
          bb=set.bounding_box();
        } 
        catch(Geometry::UnboundedSet& e) {
            throw e;
          }
        Geometry::PartitionScheme<R> ps(bb);
        int depth=16;
        Geometry::PartitionTreeSet<R> pts=Geometry::outer_approximation(set,ps,depth);
        return txt << pts;
      }
    }


    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::FiniteGrid<R>& fg)
    {
      Geometry::GridCellListSet<R> gcls(fg.grid());
      gcls.adjoin(Geometry::GridBlock<R>(fg.grid(),fg.lattice_block()));
      txt << gcls;
      return txt;
    }


    template<class R> inline
    textstream&
    operator<<(textstream& txt, const Geometry::PartitionTree<R>& pt)
    {
      for(typename Geometry::PartitionTree<R>::const_iterator iter = pt.begin(); iter!=pt.end(); ++iter) {
        txt << Geometry::Rectangle<R>(*iter);
      }
      return txt;
    }
    

    template<class BS> inline
    std::string
    summary(const Geometry::ListSet<BS>& ls) {
      std::stringstream ss;
      ss << "ListSet( size=" << ls.size() << " )";
      return ss.str();
    }

    template<class R> inline
    std::string
    summary(const Geometry::GridCellListSet<R>& gcls) {
      std::stringstream ss;
      ss << "GridCellListSet( grid=" << gcls.grid() << ", size=" << gcls.size() << " )";
      return ss.str();
    }

    template<class R> inline
    std::string
    summary(const Geometry::GridMaskSet<R>& gms) {
      std::stringstream ss;
      ss << "GridMaskSet( grid=" << gms.grid() << ", extent=" << gms.extent() << ", block=" << gms.block()
         << ", size=" << gms.size() << ", capacity=" << gms.capacity() << " )";
      return ss.str();
    }

    template<class R> inline
    std::string
    summary(const Geometry::PartitionTreeSet<R>& pts) {
      std::stringstream ss;
      ss << "PartitionTreeSet( unit_box" << pts.unit_box() << ", subdivisions=" << pts.subdivisions()
         << ", size=" << pts.size() << ", capacity=" << pts.capacity() << " )";
      return ss.str();
    }

  }
}

#endif /* ARIADNE_TEXTSTREAM_H */
