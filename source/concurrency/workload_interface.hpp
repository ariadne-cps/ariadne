/***************************************************************************
 *            concurrency/workload_interface.hpp
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

/*! \file concurrency/workload_interface.hpp
 *  \brief Interface for a workload, i.e., multiple elements to process
 */

#ifndef ARIADNE_WORKLOAD_INTERFACE_HPP
#define ARIADNE_WORKLOAD_INTERFACE_HPP

#include <functional>
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "concurrency_typedefs.hpp"

namespace Ariadne {

//! \brief Interface for a workload expressed as a stack of elements to work on, supplied with a function to process them
//! \details E: stack element type
//!          AS: optional input arguments for processing the elements; if used as output, their synchronisation
//!              in the concurrent case is up to the designer
//!          The workload handles the non-concurrent case separately, in order to unroll the tasks breadth-first: if
//!          tasks were instead enqueued to the TaskManager, a depth-first execution would be performed.
template<class E, class... AS>
class WorkloadInterface {
public:

    //! \brief Process the given elements until completion
    virtual Void process() = 0;

    //! \brief The size of the workload, i.e., the number of tasks to process
    virtual SizeType size() const = 0;

    //! \brief Append one element to process
    virtual WorkloadInterface& append(E const &e) = 0;

    //! \brief Append a list of elements to process
    virtual WorkloadInterface& append(List<E> const &es) = 0;

};

}

#endif // ARIADNE_WORKLOAD_INTERFACE_HPP