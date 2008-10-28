/***************************************************************************
 *            point.cc
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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
 
#include <sstream>
#include <string>
#include <vector>

#include "point.h"
#include "stlio.h"

namespace Ariadne {

bool 
operator<(const Point& pt1, const Point& pt2) 
{
  ARIADNE_ASSERT(pt1.dimension()==pt2.dimension());
  for(uint i=0; i!=pt1.dimension(); ++i) {
    if(pt1[i]<pt2[i]) {
      return true; 
    } else if(pt1[i]>pt2[i]) {
      return false;
    }
  }
  return false;
}
Point make_point(const std::string& str)
{
  std::vector<Float> vec;
  std::stringstream ss(str);
  read_sequence(ss,vec,'(',')',',');
  return Point(vec.size(),&vec[0]);
}

} //namespace Ariadne
