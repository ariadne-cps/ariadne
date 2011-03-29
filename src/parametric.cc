/***************************************************************************
 *            parametric.cc
 *
 *  Copyright 2010  Luca Geretti
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

#include "parametric.h"

namespace Ariadne {



ParametricOutcome::ParametricOutcome(const RealConstantSet params, const tribool value)
{
	for (RealConstantSet::const_iterator const_it = params.begin(); const_it != params.end(); ++const_it)
		_params.insert(*const_it);

	_value = value;
}

ParametricOutcome::ParametricOutcome(const ParametricOutcome& other)
{
	for (RealConstantSet::const_iterator const_it = other.getParams().begin(); const_it != other.getParams().end(); ++const_it)
		_params.insert(*const_it);

	_value = other.getOutcome();
}

ParametricOutcome&
ParametricOutcome::operator=(const ParametricOutcome& other)
{
	for (RealConstantSet::const_iterator const_it = other.getParams().begin(); const_it != other.getParams().end(); ++const_it)
		_params.insert(*const_it);

	_value = other.getOutcome();

	return *this;
}

const RealConstantSet&
ParametricOutcome::getParams() const
{
	return _params;
}

const tribool&
ParametricOutcome::getOutcome() const
{
	return _value;
}

std::ostream&
ParametricOutcome::write(std::ostream& os) const
{
	os << "(" << _params << "->" << pretty_print(_value) << ")";
	return os;
}


void
draw(std::string basename, const std::list<ParametricOutcome>& outcomes)
{
	RealConstantSet _params = outcomes.back().getParams();

	ARIADNE_ASSERT_MSG(outcomes.size() > 0, "The outcomes list is empty.");
	ARIADNE_ASSERT_MSG(_params.size() > 1, "At least two parameters are required for drawing.");

	// Plots for each couple of parameters
	for (RealConstantSet::const_iterator xparam_it = _params.begin(); xparam_it != _params.end(); ++xparam_it) {

		RealConstantSet::const_iterator yparam_it = xparam_it;
		for (++yparam_it; yparam_it != _params.end(); ++yparam_it)
		{
			std::string currentname = basename + "[" + xparam_it->name() + ","
													 + yparam_it->name() + "]";

			TextPlot trueTxt((currentname + ".true.dump").c_str());
			TextPlot falseTxt((currentname + ".false.dump").c_str());
			TextPlot indeterminateTxt((currentname + ".indeterminate.dump").c_str());

			// Sets up the figure
			Figure fig;
			Box graphics_box(2);
			graphics_box[0] = outcomes.begin()->getParams().find(*xparam_it)->value();
			graphics_box[1] = outcomes.begin()->getParams().find(*yparam_it)->value();
			array<uint> xy(2,0,1);
			fig.set_projection_map(ProjectionFunction(xy,2));

			// Adds each outcome with a dedicated fill colour
			for (std::list<ParametricOutcome>::const_iterator outcome_it = outcomes.begin();
																		  outcome_it != outcomes.end();
																		  ++outcome_it) {
				Box outcome_box(2);
				outcome_box[0] = outcome_it->getParams().find(*xparam_it)->value();
				outcome_box[1] = outcome_it->getParams().find(*yparam_it)->value();

				graphics_box[0].set_lower(min(graphics_box[0].lower(),outcome_box[0].lower()));
				graphics_box[0].set_upper(max(graphics_box[0].upper(),outcome_box[0].upper()));
				graphics_box[1].set_lower(min(graphics_box[1].lower(),outcome_box[1].lower()));
				graphics_box[1].set_upper(max(graphics_box[1].upper(),outcome_box[1].upper()));

				// Chooses the fill color and dumps the box
				tribool outcome = outcome_it->getOutcome();
				if (definitely(outcome)) {
					trueTxt.draw(outcome_box);
					fig.set_fill_colour(Colour(0.0,0.83,0.0));
				}
				else if (indeterminate(outcome)) {
					indeterminateTxt.draw(outcome_box);
					fig.set_fill_colour(Colour(1.0,1.0,0.0));
				}
				else {
					falseTxt.draw(outcome_box);
					fig.set_fill_colour(Colour(1.0,0.34,0.34));
				}

				fig.draw(outcome_box);
			}

			fig.set_bounding_box(graphics_box);
			fig.write(currentname.c_str());

			trueTxt.close();
			indeterminateTxt.close();
			falseTxt.close();
		}
	}
}

Parametric2DBisectionResults::Parametric2DBisectionResults(const std::string& filename,
							const Interval& xBounds, const Interval& yBounds,
							const unsigned& numPointsPerAxis)
{
	// Assigns
	this->filename = filename;
	this->xBounds = xBounds;
	this->yBounds = yBounds;
	this->numPointsPerAxis = numPointsPerAxis;
}


void Parametric2DBisectionResults::dump() const
{
	// Opens a file for writing
	fstream outStream;
	std::string extfilename = filename + ".dump";
	outStream.open(extfilename.c_str(), ios::out | ios::trunc);

	// For each element in the xResults
	outStream << "Varying X:\n";
	for (std::vector<std::pair<Interval,Interval> >::const_iterator it = xResults.begin(); it != xResults.end(); it++)
		outStream << it->first << it->second << ";";
	outStream << "\n";
	// For each element in the yResults
	outStream << "Varying Y:\n";
	for (std::vector<std::pair<Interval,Interval> >::const_iterator it = yResults.begin(); it != yResults.end(); it++)
		outStream << it->first << it->second << ";";

	// Closes the file
	outStream.close();
}


void Parametric2DBisectionResults::insertXValue(const std::pair<Interval,Interval>& result)
{
	// Checks that the maximum size has not been reached
	ARIADNE_ASSERT(this->xResults.size() < numPointsPerAxis);
	// Adds the value
	this->xResults.push_back(result);
}


void Parametric2DBisectionResults::insertYValue(const std::pair<Interval,Interval>& result)
{
	// Checks that the maximum size has not been reached
	ARIADNE_ASSERT(this->yResults.size() < numPointsPerAxis);
	// Adds the value
	this->yResults.push_back(result);
}


void Parametric2DBisectionResults::draw() const
{
	ARIADNE_ASSERT(this->areAllResultsAvailable());

	// Initializes the figure
	Figure fig;
	array<uint> xy(2,0,1);
	fig.set_projection_map(ProjectionFunction(xy,2));
	Box bounding_box(2,xBounds.lower(),xBounds.upper(),yBounds.lower(),yBounds.upper());
	fig.set_bounding_box(bounding_box);

	TextPlot trueTxt((filename + "-bis.true.dump").c_str());
	TextPlot falseTxt((filename + "-bis.false.dump").c_str());
	TextPlot indeterminateTxt((filename + "-bis.indeterminate.dump").c_str());
	indeterminateTxt.draw(bounding_box);
	indeterminateTxt.close();

	// Draws the bounding box in yellow
	fig.set_fill_colour(Colour(1.0,1.0,0.0));
	fig.draw(bounding_box);

	// Gets the widths and increments
	double xWidth = xBounds.width();
	double yWidth = yBounds.width();
	double xIncrement = xWidth/(numPointsPerAxis-1);
	double yIncrement = yWidth/(numPointsPerAxis-1);

	// Variables to hold the various verified and falsified intervals
	Interval leftVerified, leftFalsified, rightVerified, rightFalsified;
	Interval bottomVerified, bottomFalsified, topVerified, topFalsified;
	// The x and y spans of an identified box
	Interval boxXSpan, boxYSpan;

	// For each X result up to the last EXCLUDED
	for (unsigned i=0; i<this->xResults.size()-1; i++)
	{
		// Gets the left value for X
		double leftX = i*xIncrement+xBounds.lower();
		// Gets the right value for X
		double rightX = (i+1)*xIncrement+xBounds.lower();

		// Gets the X span of the box
		boxXSpan = Interval(leftX,rightX);

		// Gets the left and right result and saves the verified/falsified intervals
		make_lpair<Interval,Interval>(leftVerified,leftFalsified) = xResults[i];
		make_lpair<Interval,Interval>(rightVerified,rightFalsified) = xResults[i+1];

		// Draws the verified box if the intervals intersect
		if (!intersection(leftVerified,rightVerified).empty())
		{
			// Gets the Y span as ranging from the maximum of the lower values to
			// the minimum of the upper values
			boxYSpan = Interval(max(leftVerified.lower(),rightVerified.lower()),
								min(leftVerified.upper(),rightVerified.upper()));

			// Appends the set, with the desired fill color
			fig.set_fill_colour(Colour(0.0,0.83,0.0));

			Box trueBox(2,boxXSpan.lower(),boxXSpan.upper(),boxYSpan.lower(),boxYSpan.upper());

			fig.draw(trueBox);
			trueTxt.draw(trueBox);
		}

		// Draws the falsified box if the intervals intersect
		if (!intersection(leftFalsified,rightFalsified).empty())
		{

			// Gets the Y span as ranging from the maximum of the lower values to
			// the minimum of the upper values
			boxYSpan = Interval(max(leftFalsified.lower(),rightFalsified.lower()),
								min(leftFalsified.upper(),rightFalsified.upper()));

			// Appends the set, with the desired fill color
			fig.set_fill_colour(Colour(1.0,0.34,0.34));

			Box falseBox(2,boxXSpan.lower(),boxXSpan.upper(),boxYSpan.lower(),boxYSpan.upper());

			fig.draw(falseBox);
			falseTxt.draw(falseBox);
		}
	}

	// For each Y result up to the last EXCLUDED
	for (unsigned i=0; i<this->yResults.size()-1; i++)
	{
		// Gets the left value for Y
		double leftY = i*yIncrement+yBounds.lower();
		// Gets the right value for Y
		double rightY = (i+1)*yIncrement+yBounds.lower();

		// Gets the Y span of the box
		boxYSpan = Interval(leftY,rightY);

		// Gets the left and right result and saves the verified/falsified intervals
		make_lpair<Interval,Interval>(leftVerified,leftFalsified) = yResults[i];
		make_lpair<Interval,Interval>(rightVerified,rightFalsified) = yResults[i+1];

		// Draws the verified box if the intervals intersect
		if (!intersection(leftVerified,rightVerified).empty())
		{
			// Gets the X span as ranging from the maximum of the lower values to
			// the minimum of the upper values
			boxXSpan = Interval(max(leftVerified.lower(),rightVerified.lower()),
								min(leftVerified.upper(),rightVerified.upper()));

			// Appends the set, with the desired fill color
			fig.set_fill_colour(Colour(0.0,0.83,0.0));

			Box trueBox(2,boxXSpan.lower(),boxXSpan.upper(),boxYSpan.lower(),boxYSpan.upper());

			fig.draw(trueBox);
			trueTxt.draw(trueBox);
		}

		// Draws the falsified box if the intervals intersect
		if (!intersection(leftFalsified,rightFalsified).empty())
		{
			// Gets the X span as ranging from the maximum of the lower values to
			// the minimum of the upper values
			boxXSpan = Interval(max(leftFalsified.lower(),rightFalsified.lower()),
								min(leftFalsified.upper(),rightFalsified.upper()));

			// Appends the set, with the desired fill color
			fig.set_fill_colour(Colour(1.0,0.34,0.34));

			Box falseBox(2,boxXSpan.lower(),boxXSpan.upper(),boxYSpan.lower(),boxYSpan.upper());

			fig.draw(falseBox);
			falseTxt.draw(falseBox);
		}
	}

	// Draws
	fig.write(filename.c_str());
	trueTxt.close();
	falseTxt.close();
}


bool Parametric2DBisectionResults::areAllResultsAvailable() const
{
	return ((this->xResults.size() == numPointsPerAxis) && (this->yResults.size() == numPointsPerAxis));
}

} //namespace Ariadne
