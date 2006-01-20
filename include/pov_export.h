/***************************************************************************
 *            povexp.h
 *
 *  Thu Feb 10 10:33:32 2005
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
 
#ifndef _POVEXP_H
#define _POVEXP_H

#include <string>
#include <fstream>

#define POV_BBOX 40
#define POV_ACCURACY "1e-10"
#define POV_GRADIENT "70"
#define POV_COLOR "1.0, 0.8, 0.8"

void open_and_copy_pov_header(std::string &f_name, std::ofstream &os) {
  f_name+=".pov";
  os.open(f_name.c_str() , std::ios::out);

  os << "#include \"colors.inc\"" << std::endl << std::endl
     << "camera {" << std::endl << "location < " << POV_BBOX -9 
     << " , "<< POV_BBOX -9 << " , " << POV_BBOX -9 <<" >" << std::endl
     << "look_at  <0, 0, 0>"<< std::endl << "}" << std::endl
     << " light_source { <45, 40, 25> color White} "<< std::endl
   //<< " light_source { <0, 45, 0> color White} "<< std::endl
     << std::endl;
}

void close_pov_stream(std::ofstream &os) {
  os.close();
}

#endif
