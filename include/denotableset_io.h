/***************************************************************************
 *            denotableset_io.h
 *
 *  Wed Feb 16 19:32:34 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#include "rectangle.h"

#ifndef _DENOTABLESET_IO_H
#define _DENOTABLESET_IO_H

namespace Ariadne {
namespace Geometry {
	
template <class BS> class DenotableSet;

template <typename BS>
std::ostream& operator<<(std::ostream &os, 
				const DenotableSet<BS> &A){

	if (A.size() >0 )
		os << "BasicSet[0]=" << A[0];

	for (size_t i=1; i<A.size(); i++) 
			os << std::endl << "BasicSet["<< i <<"]=" << A[i];
	

	return os;	
}

template <typename BS>
std::istream& operator>>(std::istream &is, 
				DenotableSet<BS> &A)
{
    std::vector<BS>& vec(A._vector);
    is >> vec;
    
    if(vec.size()==0) {
	A._dim = 0;
    }
    else {
	A._dim=vec[0].dim();
    }

    return is;
}


template <typename BSE>
class DenotableSetExporter {
	
	public:
		typedef BSE Exporter;
		typedef typename Exporter::BasicSet BasicSet;
		typedef typename BasicSet::State State;
		typedef typename BasicSet::Real Real;

                typedef DenotableSet<BasicSet> DenotableSet;
		typedef Rectangle<State> Rectangle;
	
	private:
		
		Exporter _e;
	
	public:
	
		DenotableSetExporter() {}
			
		DenotableSetExporter(const Exporter &e): _e(e) {}
			
		DenotableSetExporter(const std::string &file_name, 
			const Rectangle &bbox): _e(file_name,bbox) {}
			
		DenotableSetExporter(const std::string &file_name, 
				const Rectangle &bbox, std::vector<uint> dims): 
				_e(file_name, bbox, dims) {}
			 
		DenotableSetExporter(std::ofstream os, Rectangle bbox): _e(os,bbox) {}
			
		DenotableSetExporter(const std::string &file_name): _e(bbox) {}
			
		DenotableSetExporter(const std::string &file_name,
					std::vector<uint> dims): 
				_e(file_name, dims) {}
					
		DenotableSetExporter(std::ofstream os): _e(os) {}
			
		void export_denotableset(const DenotableSet &A) {
				
			for (size_t i=0; i< A.size(); i++) {
				(this->_e).export_basicset(A[i]); 
			}
		}
					
		void export_denotableset(const DenotableSet &A, std::vector<uint> dims) {
				
			for (size_t i=0; i< A.size(); i++) {
				(this->_e).export_basicset(A[i], dims); 
			}
		}
};	

}
}

#endif /* _DENOTABLESET_IO_H */
