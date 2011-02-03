/***************************************************************************
 *            parametric.h
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

/*! \file parametric2d.h
 *  \brief Data structures for handling parametric analysis in 2D.
 */

#include "ariadne.h"

#ifndef ARIADNE_PARAMETRIC_H
#define ARIADNE_PARAMETRIC_H

namespace Ariadne {

typedef std::set<RealConstant,ConstantNameComparator<Real> > RealConstantSet;

/**
 * \brief The data structure for the verification outcome over a configuration of parameters (i.e. constants of a system)
 */
struct ParametricVerificationOutcome
{
private:

	RealConstantSet _params;

	tribool _value;

public:

	ParametricVerificationOutcome(const RealConstantSet params, const tribool value);
	ParametricVerificationOutcome(const ParametricVerificationOutcome& other);

	ParametricVerificationOutcome& operator=(const ParametricVerificationOutcome& other);

	virtual std::ostream& write(std::ostream&) const;

	/** The parameters */
	const RealConstantSet& getParams() const;
	/** The value of the verification */
	const tribool& getValue() const;

};

inline std::ostream& operator<<(std::ostream& os, const ParametricVerificationOutcome& outcome) {
    return outcome.write(os); }

/**
 * \brief The list of parametric verification outcomes
 */
struct ParametricVerificationOutcomeList
{
private:

	RealConstantSet _params;

	std::list<ParametricVerificationOutcome> _outcomes;

public:

	ParametricVerificationOutcomeList(const RealConstantSet& params);

	virtual std::ostream& write(std::ostream&) const;

	/** Inserts a result */
	void push_back(const ParametricVerificationOutcome& result);

	/** The current outcomes */
	const std::list<ParametricVerificationOutcome>& getOutcomes() const;

};

inline std::ostream& operator<<(std::ostream& os, const ParametricVerificationOutcomeList& outcomeList) {
    return outcomeList.write(os); }

/**
 * \brief A data structure for the results on a parametric 2D analysis using bisection.
 */
struct Parametric2DBisectionResults
{
private:

	/**
	 * The file name used for dumping and plotting.
	 */
	std::string filename;

	/**
	 * The bounds for the X coordinate.
	 */
	Interval xBounds;

	/**
	 * The bounds for the Y coordinate.
	 */
	Interval yBounds;

	/**
	 * The number of points per axis used for the analysis.
	 */
	unsigned numPointsPerAxis;

	/**
	 * The vector of the results for fixed X values.
	 */
	std::vector<std::pair<Interval,Interval> > xResults;

	/**
	 * The vector of the results for fixed Y values.
	 */
	std::vector<std::pair<Interval,Interval> > yResults;

public:

	/**
	 * \brief Constructs from parameters.
	 *
	 * @param filename The file name to use for dumping and plotting.
	 * @param xBounds The X bounds.
	 * @param yBounds The Y bounds.
	 * @param numPointsPerAxis The number of points per axis used.
	 */
	Parametric2DBisectionResults(const std::string& filename,
								const Interval& xBounds, const Interval& yBounds,
								const unsigned& numPointsPerAxis);

	/**
	 * \brief Dumps the data of the results into a file.
	 */
	void dump() const;

	/**
	 * \brief Adds a value for the analysis on X.
	 *
	 * @param result The pair of the value + the verified-falsified intervals.
	 */
	void insertXValue(const std::pair<Interval,Interval>& result);

	/**
	 * \brief Adds a value for the analysis on Y.
	 *
	 * @param result The pair of the value + the verified-falsified intervals.
	 */
	void insertYValue(const std::pair<Interval,Interval>& result);

	/**
	 * \brief Draws the result in the current folder.
	 */
	void draw() const;

	/**
	 * \brief Checks if all the results are available.
	 *
	 * @return True if the sizes of the results match with the maximum ones, false otherwise.
	 */
	bool areAllResultsAvailable() const;

};

} // namespace Ariadne

#endif // ARIADNE_PARAMETRIC_H
