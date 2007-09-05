/****************************************************************************
 *            chompfstream.cc
 *
 *  Copyright  2007  Pieter Collins, Davide Bresolin
 *  Pieter.Collins@cwi.nl, bresolin@sci.univr.it
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

#include "output/chompstream.h"

namespace Ariadne { 



Output::chompstream::chompstream(std::ostream& os)
 : _os_ptr(&os)
{
}


void
Output::chompstream::redirect(std::ostream& os)
{
  this->_os_ptr=&os;
}


Output::chompstream::~chompstream() {
}

      

Output::chompstream& 
Output::operator<<(chompstream& chomps, const char* str) 
{
  std::ostream& os=*chomps._os_ptr;
  os << str; 
  return chomps;
}

    
Output::chompstream& 
Output::operator<<(chompstream& chomps, const Combinatoric::LatticeCell& lc) 
{
  std::ostream& os=*chomps._os_ptr;
  if(lc.dimension()>0) {
    os << "(" << lc.lower_bound(0);
    for(dimension_type i=1; i!=lc.dimension(); ++i) {
      os << "," << lc.lower_bound(i);
    }
    os << ")";
  } else {
    os << "()";
  }
  return chomps;
}


Output::chompstream& 
Output::operator<<(chompstream& chomps, const Combinatoric::LatticeMaskSet& lms) 
{
  for(Combinatoric::LatticeMaskSet::const_iterator cell_iter=lms.begin(); cell_iter!=lms.end(); ++cell_iter) {
    chomps << *cell_iter << "\n";
  }
  return chomps;
}


Output::chompstream& 
Output::operator<<(chompstream& chomps, const Combinatoric::LatticeMultiMap& lmm) 
{
  Combinatoric::LatticeCell lc(lmm.argument_dimension());
  Combinatoric::LatticeCellListSet lcls(lmm.result_dimension());
  for(Combinatoric::LatticeMultiMap::const_iterator iter=lmm.begin(); iter!=lmm.end(); ++iter) {
    lc=iter->first;
    lcls=iter->second;
    chomps << lc << " -> { ";
    for(Combinatoric::LatticeCellListSet::const_iterator lcls_iter=lcls.begin(); 
        lcls_iter!=lcls.end(); ++lcls_iter) 
    {
      chomps << *lcls_iter << " ";
    }
    chomps << "}\n";
  }
  return chomps;
}


    
}
