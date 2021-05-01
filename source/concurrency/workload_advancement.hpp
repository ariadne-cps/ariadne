/***************************************************************************
 *            concurrency/workload_advancement.hpp
 *
 *  Copyright  2007-21  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file concurrency/workload_advancement.hpp
 *  \brief Synchronised class to manage the status of multiple elements to process
 */

#ifndef ARIADNE_WORKLOAD_ADVANCEMENT_HPP
#define ARIADNE_WORKLOAD_ADVANCEMENT_HPP

#include <algorithm>
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "concurrency_typedefs.hpp"

namespace Ariadne {

//! \brief Synchronised class to manage the status of multiple elements to process
class WorkloadAdvancement {
  public:
    WorkloadAdvancement(SizeType initial = 0);

    //! \brief The elements waiting to be processed
    SizeType waiting() const;
    //! \brief The elements under processing
    SizeType processing() const;
    //! \brief The completed elements
    SizeType completed() const;

    //! \brief Add n elements to waiting
    Void add_to_waiting(SizeType n = 1);
    //! \brief Move n waiting to processing
    Void add_to_processing(SizeType n = 1);
    //! \brief Move n processing to completed
    Void add_to_completed(SizeType n = 1);

    //! \brief The rate of completion r (0<=r<=1) related to the progress
    Double completion_rate() const;

    //! \brief If no other processing remains
    Bool has_finished() const;

  private:
    SizeType _num_waiting;
    SizeType _num_processing;
    SizeType _num_completed;

    Mutex mutable _mux;
};

} // namespace Ariadne

#endif // ARIADNE_WORKLOAD_ADVANCEMENT_HPP
