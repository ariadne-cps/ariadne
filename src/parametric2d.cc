/***************************************************************************
 *            parametric2d.cc
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

#include "parametric2d.h"

namespace Ariadne {


Parametric2DAnalysisResults::Parametric2DAnalysisResults(const std::string& filename,
							const Interval& xBounds, const Interval& yBounds,
							const unsigned& numPointsPerAxis)
{
	// Assigns
	this->filename = filename;
	this->xBounds = xBounds;
	this->yBounds = yBounds;
	this->numPointsPerAxis = numPointsPerAxis;
}


void Parametric2DAnalysisResults::dump() const
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


void Parametric2DAnalysisResults::insertXValue(const std::pair<Interval,Interval>& result)
{
	// Checks that the maximum size has not been reached
	assert(this->xResults.size() < numPointsPerAxis);
	// Adds the value
	this->xResults.push_back(result);
}


void Parametric2DAnalysisResults::insertYValue(const std::pair<Interval,Interval>& result)
{
	// Checks that the maximum size has not been reached
	assert(this->yResults.size() < numPointsPerAxis);
	// Adds the value
	this->yResults.push_back(result);
}


void Parametric2DAnalysisResults::draw() const
{
	assert(this->areAllResultsAvailable());

	// Initializes the figure
	Figure fig;
	array<uint> xy(2,0,1);
	fig.set_projection_map(ProjectionFunction(xy,2));
	Box bounding_box(2,xBounds.lower(),xBounds.upper(),yBounds.lower(),yBounds.upper());
	fig.set_bounding_box(bounding_box);

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

		// Draws the verified box (if neither the left or right intervals are empty)
		if (!leftVerified.empty() && !rightVerified.empty())
		{
			// Gets the Y span as ranging from the maximum of the lower values to
			// the minimum of the upper values
			boxYSpan = Interval(max(leftVerified.lower(),rightVerified.lower()),
								min(leftVerified.upper(),rightVerified.upper()));

			// Appends the set, with the desired fill color
			fig.set_fill_colour(Colour(0.0,0.83,0.0));
			fig.draw(Box(2,boxXSpan.lower(),boxXSpan.upper(),boxYSpan.lower(),boxYSpan.upper()));
		}

		// Draws the falsified box (if neither the left or right intervals are empty)
		if (!leftFalsified.empty() && !rightFalsified.empty())
		{
			// Gets the Y span as ranging from the maximum of the lower values to
			// the minimum of the upper values
			boxYSpan = Interval(max(leftFalsified.lower(),rightFalsified.lower()),
								min(leftFalsified.upper(),rightFalsified.upper()));

			// Appends the set, with the desired fill color
			fig.set_fill_colour(Colour(1.0,0.34,0.34));
			fig.draw(Box(2,boxXSpan.lower(),boxXSpan.upper(),boxYSpan.lower(),boxYSpan.upper()));
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

		// Draws the verified box (if neither the left or right intervals are empty)
		if (!leftVerified.empty() && !rightVerified.empty())
		{
			// Gets the X span as ranging from the maximum of the lower values to
			// the minimum of the upper values
			boxXSpan = Interval(max(leftVerified.lower(),rightVerified.lower()),
								min(leftVerified.upper(),rightVerified.upper()));

			// Appends the set, with the desired fill color
			fig.set_fill_colour(Colour(0.0,0.83,0.0));
			fig.draw(Box(2,boxXSpan.lower(),boxXSpan.upper(),boxYSpan.lower(),boxYSpan.upper()));
		}

		// Draws the falsified box (if neither the left or right intervals are empty)
		if (!leftFalsified.empty() && !rightFalsified.empty())
		{
			// Gets the X span as ranging from the maximum of the lower values to
			// the minimum of the upper values
			boxXSpan = Interval(max(leftFalsified.lower(),rightFalsified.lower()),
								min(leftFalsified.upper(),rightFalsified.upper()));

			// Appends the set, with the desired fill color
			fig.set_fill_colour(Colour(1.0,0.34,0.34));
			fig.draw(Box(2,boxXSpan.lower(),boxXSpan.upper(),boxYSpan.lower(),boxYSpan.upper()));
		}
	}

	// Draws
	fig.write(filename.c_str());
}


bool Parametric2DAnalysisResults::areAllResultsAvailable() const
{
	return ((this->xResults.size() == numPointsPerAxis) && (this->yResults.size() == numPointsPerAxis));
}

} //namespace Ariadne
