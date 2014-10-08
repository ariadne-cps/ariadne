/***************************************************************************
 *            textplot.cc
 *
 *  Copyright 2009  Davide Bresolin
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

#include "functional.h"
#include "config.h"

#include "macros.h"
#include "stlio.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "point.h"
#include "box.h"
#include "curve.h"
//#include "polytope.h"
//#include "zonotope.h"
#include "textplot.h"
#include "grid_set.h"

namespace Ariadne {


TextPlot::~TextPlot()
{
    this->_fstream.close();
}


TextPlot::TextPlot()
    : _fstream()
{
}

TextPlot::TextPlot(const char* cfilename)
{
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), ios::out | ios::trunc);
}


TextPlot::TextPlot(const char* cfilename, ios_base::openmode mode)
{
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), mode);
}


void TextPlot::open(const char* cfilename)
{
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), ios::out | ios::trunc);
}


void TextPlot::open(const char* cfilename, ios_base::openmode mode)
{
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), mode);
}


inline std::ostream& operator<<(std::ostream& os, const DrawableInterface& sh) { return sh.write(os); }

void TextPlot::draw(const DrawableInterface& shape) {
    if(dynamic_cast<const ExactPoint*>(&shape)) {
        this->draw(dynamic_cast<const ExactPoint&>(shape));
    } else if(dynamic_cast<const ExactBox*>(&shape)) {
        this->draw(dynamic_cast<const ExactBox&>(shape));
//    } else if(dynamic_cast<const Polytope*>(&shape)) {
//        this->draw(dynamic_cast<const Polytope&>(shape));
    } else if(dynamic_cast<const InterpolatedCurve*>(&shape)) {
        this->draw(dynamic_cast<const InterpolatedCurve&>(shape));
    } else if(dynamic_cast<const GridTreeSubset*>(&shape)) {
        this->draw(dynamic_cast<const GridTreeSubset&>(shape));
    } else {
        ARIADNE_THROW(std::runtime_error,"TextPlot::draw(const DrawableInterface&)","Unrecognised shape "<<shape<<" "<<typeid(shape).name());
    }
}

void TextPlot::_draw(const std::vector<ExactPoint>& pts) {
    for(std::vector<ExactPoint>::const_iterator iter = pts.begin() ; iter != pts.end() ; iter++) {
        this->draw(*iter);
    }
    this->_fstream << std::endl;
}

void TextPlot::draw(const ExactPoint& pt) {
    for(uint i = 0; i < pt.dimension(); i++) {
        this->_fstream << approx_cast<double>(pt[i]) << " ";
    }
    this->_fstream << std::endl;
}


void TextPlot::draw(const ExactBox& bx) {
    this->_draw(bx.vertices());
}

//void TextPlot::draw(const Polytope& p) {
//    this->_draw(p.vertices());
//}

void TextPlot::draw(const InterpolatedCurve& c) {
    for(InterpolatedCurve::const_iterator iter = c.begin() ; iter != c.end() ; ++iter) {
        this->draw(ApproximatePoint(iter->second));
    }
    this->_fstream << std::endl;
}

void TextPlot::draw(const GridTreeSubset& gts) {
    for(GridTreeSubset::const_iterator iter = gts.begin() ; iter != gts.end() ; ++iter) {
        this->draw(iter->box());
    }
    this->_fstream << std::endl;
}

void TextPlot::close() {
    this->_fstream.close();
}


} // Namespace Ariadne


