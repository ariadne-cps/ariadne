/***************************************************************************
 *            concurrency/buffer.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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
#include "numeric/logical.hpp"

namespace Ariadne {

template<class E> class Buffer
{
  public:
    Buffer(unsigned int size) : _size(size), _end_consuming(false) { }

    void push(E const& e) {
        std::unique_lock<std::mutex> locker(mu);
        cond.wait(locker, [this](){return _buffer.size() < _size;});
        _buffer.push(e);
        locker.unlock();
        cond.notify_all();
    }
    E pop() {
        std::unique_lock<std::mutex> locker(mu);
        cond.wait(locker, [this](){return _buffer.size() > 0 || _end_consuming;});

        if (_end_consuming) throw std::exception();
        E back = _buffer.front();
        _buffer.pop();
        locker.unlock();
        cond.notify_all();
        return back;
    }

    void end_consuming() {
        _end_consuming = true;
        cond.notify_all();
    }

    bool ended_consuming() const {
        return _end_consuming;
    }

private:
    std::mutex mu;
    std::condition_variable cond;
    std::queue<E> _buffer;
    const unsigned int _size;
    std::atomic<bool> _end_consuming;
};

} // namespace Ariadne

#endif // ARIADNE_BUFFER_HPP
