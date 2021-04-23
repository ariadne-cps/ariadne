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
#include "numeric/logical.hpp"

namespace Ariadne {

//! \brief Exception useful when the buffer is allowed to stay in the receiving condition
class BufferStoppedConsumingException : public std::exception { };

template<class E> class Buffer
{
  public:
    Buffer(SizeType capacity) : _capacity(capacity), _stop_consuming(false) { ARIADNE_PRECONDITION(capacity>0); }

    //! \brief Push an object into the buffer
    //! \details Will block if the capacity has been reached
    void push(E const& e) {
        std::unique_lock<std::mutex> locker(mux);
        cond.wait(locker, [this](){return _queue.size() < _capacity;});
        _queue.push(e);
        locker.unlock();
        cond.notify_all();
    }

    //! \brief Pulls an object from the buffer
    //! \details Will block if the capacity is zero
    E pull() {
        std::unique_lock<std::mutex> locker(mux);
        cond.wait(locker, [this](){return _queue.size() > 0 || _stop_consuming;});
        if (_stop_consuming) { throw BufferStoppedConsumingException(); }
        E back = _queue.front();
        _queue.pop();
        locker.unlock();
        cond.notify_all();
        return back;
    }

    //! \brief The current size of the queue
    SizeType size() const {
        std::lock_guard<std::mutex> locker(mux);
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

    //! \brief Trigger a stop to consuming in the case that the buffer was in the waiting state for input
    void stop_consuming() {
        _stop_consuming = true;
        cond.notify_all();
    }

private:
    mutable std::mutex mux;
    std::condition_variable cond;
    std::queue<E> _queue;
    std::atomic<SizeType> _capacity;
    Bool _stop_consuming;
};

} // namespace Ariadne

#endif // ARIADNE_BUFFER_HPP
