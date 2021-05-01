/***************************************************************************
 *            concurrency/buffer.hpp
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

/*! \file concurrency/buffer.hpp
 *  \brief A multiple-thread-safe queue usable as a buffer.
 */

#ifndef ARIADNE_BUFFER_HPP
#define ARIADNE_BUFFER_HPP

#include <utility>
#include <mutex>
#include <condition_variable>
#include <queue>
#include "utility/typedefs.hpp"
#include "utility/macros.hpp"
#include "concurrency_typedefs.hpp"

namespace Ariadne {

//! \brief Exception useful when the buffer is allowed to stay in the receiving condition
class BufferInterruptPullingException : public std::exception { };

//! \brief A class for handling a buffer
template<class E> class Buffer
{
  public:
    Buffer(SizeType capacity) : _capacity(capacity), _interrupt(false) { ARIADNE_PRECONDITION(capacity > 0); }

    //! \brief Push an object into the buffer
    //! \details Will block if the capacity has been reached
    void push(E const& e) {
        UniqueLock<Mutex> locker(mux);
        cond.wait(locker, [this](){return _queue.size() < _capacity;});
        _queue.push(e);
        cond.notify_all();
    }

    //! \brief Pulls an object from the buffer
    //! \details Will block if the capacity is zero
    E pull() {
        UniqueLock<Mutex> locker(mux);
        cond.wait(locker, [this](){return not _queue.empty() || _interrupt;});
        if (_interrupt and _queue.empty()) { _interrupt = false; throw BufferInterruptPullingException(); }
        E back = _queue.front();
        _queue.pop();
        cond.notify_all();
        return back;
    }

    //! \brief The current size of the queue
    SizeType size() const {
        LockGuard<Mutex> locker(mux);
        return _queue.size();
    }

    //! \brief The maximum size for the queue
    SizeType capacity() const {
        return _capacity;
    }

    //! \brief Change the capacity
    Void set_capacity(SizeType capacity) {
        ARIADNE_PRECONDITION(capacity>0);
        ARIADNE_ASSERT_MSG(capacity>=size(),"Reducing capacity below currenty buffer size is not allowed.");
        _capacity = capacity;
    }

    //! \brief Interrupt consuming in the case that the queue is empty and the buffer in the waiting state for input
    //! \details Needs to
    void interrupt_consuming() {
        _interrupt = true;
        cond.notify_all();
    }

private:
    mutable Mutex mux;
    ConditionVariable cond;
    std::queue<E> _queue;
    std::atomic<SizeType> _capacity;
    Bool _interrupt;
};

} // namespace Ariadne

#endif // ARIADNE_BUFFER_HPP
