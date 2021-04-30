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
//!          AS: optional input arguments for processing the elements; if used as output, their synchronisation
//!              in the concurrent case is up to the designer
template<class E, class... AS>
class Workload {
  public:
    //! \brief Access for pushing elements to the stack when processing the workload
    class StackAccess {
        friend Workload;
      protected:
        StackAccess(Workload& parent) : _load(parent) { }
      public:
        Void push(E const &e) { _load.add(e); }
      private:
        Workload& _load;
    };
  public:
    using FunctionType = std::function<Void(Workload<E,AS...>::StackAccess&, E, AS...)>;
    using BoundFunctionType = std::function<Void(E)>;

    Workload(FunctionType f, AS... as) : _sa(StackAccess(*this)), _f(std::bind(
                std::forward<FunctionType const>(f),
                std::forward<Workload<E,AS...>::StackAccess const&>(_sa),
                std::placeholders::_1,
                std::forward<AS>(as)...)) { }

    Void process() {
        ARIADNE_PRECONDITION(not _elements.empty());
        while (not _elements.empty()) {
            auto e = _elements.back();
            _elements.pop_back();
            (*e)();
        }
    }

    //! \brief Add one element to process
    Void add(E const &e) {
        LockGuard<Mutex> lock(_mux);
        _elements.push_back(std::make_shared<PackagedTask<Void()>>(std::bind(
                std::forward<BoundFunctionType const>(_f),
                std::forward<E const&>(e))));
    }

    //! \brief The size of the stack
    SizeType size() const { return _elements.size(); }

  private:
    List<SharedPointer<PackagedTask<Void()>>> _elements;
    StackAccess const _sa;
    BoundFunctionType const _f;
    Mutex _mux;
};

}

#endif // ARIADNE_WORKLOAD_HPP