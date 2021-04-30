/***************************************************************************
 *            concurrency/workload.hpp
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

/*! \file concurrency/workload.hpp
 *  \brief A stack of objects to work on, supplied with a function to process them
 */

#ifndef ARIADNE_WORKLOAD_HPP
#define ARIADNE_WORKLOAD_HPP

#include <functional>
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "concurrency_typedefs.hpp"

namespace Ariadne {

//! \brief A workload expressed as a stack of elements to work on, supplied with a function to process them
//! \details E: stack element type
//!          AS: optional input arguments for processing the elements; if passed by
//!              reference, the objects should be synchronised but this is left to the implementation
template<class E, class... AS>
class Workload {
  public:
    //! \brief Access for pushing elements to the stack when processing the workload
    class StackAccess {
        friend Workload;
      protected:
        StackAccess(Workload& parent) : _load(parent) { }
      public:
        Void push(E const &e) { LockGuard<Mutex> lock(_mux); _load._elements.push_back(e); }
      private:
        Workload& _load;
        Mutex _mux;
    };
  public:
    using FunctionType = std::function<Void(Workload<E,AS...>::StackAccess&, E, AS...)>;

    Workload(FunctionType f) : _f(f) { }

    Void process(AS... args) {
        ARIADNE_PRECONDITION(not _elements.empty());
        auto access = Workload::StackAccess(*this);
        while (not _elements.empty()) {
            auto e = _elements.back();
            _elements.pop_back();
            _f(access,e,args...);
        }
    }

    //! \brief Add one element to process
    Void add(E const &e) { _elements.push_back(e); }
    //! \brief Add a list of elements to process
    Void add(List<E> const &es) { _elements.append(es); }

    //! \brief The size of the stack
    SizeType size() const { return _elements.size(); }

  private:
    List<E> _elements;
    FunctionType const _f;
};

}

#endif // ARIADNE_WORKLOAD_HPP