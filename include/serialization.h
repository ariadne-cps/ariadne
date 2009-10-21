/****************************************************************************
 *            serialization.h
 *
 *  Copyright  2008  Pieter Collins
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

#ifndef ARIADNE_SERIALIZATION_H
#define ARIADNE_SERIALIZATION_H

/*! \file serialization.h
 *  \brief Reading and writing to a boost archive.
 */
 
//#define BOOST_NO_INT64_T

#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

//#include <boost/serialization/vector.hpp>
//#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>

namespace Ariadne {
 
using boost::archive::text_oarchive;
using boost::archive::text_iarchive;

class InputArchive : public text_iarchive { };
class OutputArchive : public text_oarchive { };

}


#endif /* ARIADNE_SERIALIZATION_H */
