/***************************************************************************
 *            matlabexp.h
 *
 *  Wed Feb 16 18:17:41 2005
 *  Copyright  2005  Alberto Casagrande, Pierpaolo Murrieri
 *  casagrande@dimi.uniud.it, murrieri@parades.rm.cnr.it
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
 
#ifndef _MATLABEXP_H
#define _MATLABEXP_H

#include <string>
#include <fstream>

void open_matlab_stream(std::string &f_name, std::ofstream &os) {
  f_name+=".m";
  os.open(f_name.c_str() , std::ios::out);
}

void close_mathlab_stream(std::ofstream &os) {
  os.close();
}

void open_OneDim_matlab_stream(std::string &f_name, std::ofstream &os) {
  f_name+=".m";
  os.open(f_name.c_str() , std::ios::out);
  os<< "data = [ ";
}

void close_OneDim_mathlab_stream(std::ofstream &os) {
  os<< "] "<<std::endl;
  os.close();
}


#endif /* _MATLABEXP_H */
