/***************************************************************************
 *            textplot.cpp
 *
 *  Copyright 2009--17  Davide Bresolin
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../utility/macros.hpp"
#include "../utility/stlio.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../function/function.hpp"
#include "../geometry/point.hpp"
#include "../geometry/box.hpp"
#include "../geometry/curve.hpp"
//#include "../geometry/polytope.hpp"
//#include "../geometry/zonotope.hpp"
#include "../output/textplot.hpp"
#include "../geometry/grid_paving.hpp"

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
    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), std::ios::out | std::ios::trunc);
}


TextPlot::TextPlot(const char* cfilename, std::ios_base::openmode mode)
{
    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), mode);
}


Void TextPlot::open(const char* cfilename)
{
    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), std::ios::out | std::ios::trunc);
}


Void TextPlot::open(const char* cfilename, std::ios_base::openmode mode)
{
    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), mode);
}


inline OutputStream& operator<<(OutputStream& os, const DrawableInterface& sh) { return sh.write(os); }

TextPlot& TextPlot::draw(const DrawableInterface& shape) {
    if(dynamic_cast<const ExactPoint*>(&shape)) {
        this->draw(dynamic_cast<const ExactPoint&>(shape));
    } else if(dynamic_cast<const ExactBoxType*>(&shape)) {
        this->draw(dynamic_cast<const ExactBoxType&>(shape));
//    } else if(dynamic_cast<const Polytope*>(&shape)) {
//        this->draw(dynamic_cast<const Polytope&>(shape));
    } else if(dynamic_cast<const InterpolatedCurve*>(&shape)) {
        this->draw(dynamic_cast<const InterpolatedCurve&>(shape));
    } else if(dynamic_cast<const GridTreeSubpaving*>(&shape)) {
        this->draw(dynamic_cast<const GridTreeSubpaving&>(shape));
    } else {
        ARIADNE_THROW(std::runtime_error,"TextPlot::draw(const DrawableInterface&)","Unrecognised shape "<<shape<<" "<<typeid(shape).name());
    }
    return *this;
}

TextPlot& TextPlot::_draw(const std::vector<ExactPoint>& pts) {
    for(std::vector<ExactPoint>::const_iterator iter = pts.begin() ; iter != pts.end() ; iter++) {
        this->draw(*iter);
    }
    this->_fstream << std::endl;
    return *this;
}

TextPlot& TextPlot::draw(const ExactPoint& pt) {
    for(Nat i = 0; i < pt.dimension(); i++) {
        this->_fstream << numeric_cast<double>(pt[i]) << " ";
    }
    this->_fstream << std::endl;
    return *this;
}


TextPlot& TextPlot::draw(const ExactBoxType& bx) {
    this->_draw(bx.vertices());
    return *this;
}

//TextPlot& TextPlot::draw(const Polytope& p) {
//    this->_draw(p.vertices());
//}

TextPlot& TextPlot::draw(const InterpolatedCurve& c) {
    for(InterpolatedCurve::ConstIterator iter = c.begin() ; iter != c.end() ; ++iter) {
        this->draw(ApproximatePoint(iter->second));
    }
    this->_fstream << std::endl;
    return *this;
}

TextPlot& TextPlot::draw(const GridTreeSubpaving& gts) {
    for(GridTreeSubpaving::ConstIterator iter = gts.begin() ; iter != gts.end() ; ++iter) {
        this->draw(iter->box());
    }
    this->_fstream << std::endl;
    return *this;
}

Void TextPlot::close() {
    this->_fstream.close();
}


} // Namespace Ariadne


