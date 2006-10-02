/****************************************************************************
 *            chompfstream.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

/*! \file chompfstream.h
 *  \brief Output to CHomP computational homology programs
 */
 
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include "../combinatoric/lattice_set.h"
#include "../combinatoric/lattice_map.h"

namespace Ariadne {
  namespace CHomP {
    /*
    class chompfstream;
    
    chompfstream& operator<<(chompfstream&, const std::string&);
    chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeMaskSet&);
    chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeCellListSet&);
    chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeCell&);
 */
    
    class chompfstream {
     public:
      chompfstream() : _ofs() { }
      chompfstream(const char* fn) : _ofs(fn) { }
      ~chompfstream() { this->close(); }

      void open(const char* fn) { _ofs.open(fn); }
      void close() { _ofs.close(); }
     private:
      friend chompfstream& operator<<(chompfstream&, const char*);
      friend chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeCell&);
      friend chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeMaskSet&);
      friend chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeMultiMap&);
     private:
      std::ofstream _ofs;
    };

    
    inline chompfstream& operator<<(chompfstream& cfs, const char* str) {
      cfs._ofs << str; return cfs;
    }
    
    inline chompfstream& operator<<(chompfstream& cfs, const Combinatoric::LatticeCell& lc) {
      std::ofstream& ofs=cfs._ofs;
      if(lc.dimension()>0) {
        ofs << "(" << lc.lower_bound(0);
        for(dimension_type i=1; i!=lc.dimension(); ++i) {
          ofs << "," << lc.lower_bound(i);
        }
        ofs << ")";
      } else {
        ofs << "()";
      }
      return cfs;
    }
    
    inline chompfstream& operator<<(chompfstream& cfs, const Combinatoric::LatticeMaskSet& lms) {
      Combinatoric::LatticeCell lc;
      for(Combinatoric::LatticeMaskSet::const_iterator iter=lms.begin(); iter!=lms.end(); ++iter) {
        lc=*iter;
        cfs << lc;
        cfs._ofs << "\n";
      }
      return cfs;
    }
    
    inline chompfstream& operator<<(chompfstream& cfs, const Combinatoric::LatticeMultiMap& lmm) {
      Combinatoric::LatticeCell lc;
      Combinatoric::LatticeCellListSet lcls;
      for(Combinatoric::LatticeMultiMap::const_iterator iter=lmm.begin(); iter!=lmm.end(); ++iter) {
        lc=iter->first;
        lcls=iter->second;
        cfs << lc << " -> { ";
        for(Combinatoric::LatticeCellListSet::const_iterator lcls_iter=lcls.begin(); 
            lcls_iter!=lcls.end(); ++lcls_iter) 
        {
          cfs << *lcls_iter << " ";
        }
        cfs << "}\n";
      }
      return cfs;
    }
    
  }
}
