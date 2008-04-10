/****************************************************************************
 *            epsstream.inline.h
 *
 *  Copyright  2005-7  Pieter Collins
 *
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

namespace Ariadne {

inline 
Output::epsstream&
Output::operator<<(epsstream& eps, const char* s) 
{
  eps.ostream() << s; 
  return eps;
}





inline Output::epsfstream::epsfstream()
  : epsstream(), 
    _ofs_ptr(new std::ofstream())  
{ 
  epsstream::redirect(*this->_ofs_ptr); 
}

template<class R> inline 
void 
Output::epsfstream::open(const char* fn, const Geometry::Box<R>& bbox) 
{
  PlanarProjectionMap p_map(bbox.dimension(),0,1);
  Rectangle2d bbox2d=p_map(bbox);
  this->open(fn,bbox2d,p_map);
}

template<class R> inline 
void 
Output::epsfstream::open(const char* fn, const Geometry::Box<R>& bbox, uint ix, uint iy) 
{
  PlanarProjectionMap p_map(bbox.dimension(),ix,iy);
  Rectangle2d bbox2d=p_map(bbox);
  this->open(fn,bbox2d,p_map);
}

template<class R> inline 
void 
Output::epsfstream::open(const char* fn, const Geometry::Box<R>& bbox, 
                         const PlanarProjectionMap& p_map)
{
  Rectangle2d bbox2d=p_map(bbox);
  this->open(fn,bbox2d,p_map);
}

inline
void 
Output::epsfstream::close() 
{ 
  this->write_trailer();
  this->_ofs_ptr->close(); 
}

inline 
Output::epsfstream::~epsfstream() 
{ 
  this->close();
  delete this->_ofs_ptr; 
}


} //namespace Ariadne
