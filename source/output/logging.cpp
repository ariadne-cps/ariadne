/***************************************************************************
 *            output/logging.cpp
 *
 *  Copyright  2007-20  Pieter Collins
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

#include <iostream>
#include <cassert>
#include <sys/ioctl.h>
#include <unistd.h>

#include "../config.hpp"
#include "../output/logging.hpp"
#include "../utility/writable.hpp"
#include "../utility/macros.hpp"

namespace Ariadne {

std::ofstream log_file_stream;

void redirect_log(const char* filename)
{
    if(log_file_stream.is_open()) {
        log_file_stream.close();
    }
    log_file_stream.open(filename);
    std::clog.rdbuf( log_file_stream.rdbuf() );
}

OutputStream& operator<<(OutputStream& os, const WritableInterface& w);

void write_error(OutputStream& os, const char* i, const char* c, const char* t) {
    os << "Error in dynamic_handle_cast: cannot convert object of static type " << c << " and dynamic type " << i << " to type " << t << "\n";
}

void write_error(OutputStream& os, const WritableInterface* w, const char* i, const char* c, const char* t) {
    os << "Error in dynamic_handle_cast:" << std::flush;
    os << " cannot convert "; assert(w); w->_write(os); os << std::flush;
    os << " of static type " << c << " and dynamic type " << i << " to type " << t << std::endl;
}

std::string very_pretty_function(std::string msg) {
    size_t ariadne_pos = std::string::npos;
    while ((ariadne_pos  = msg.find("Ariadne::")) != std::string::npos) msg.erase(ariadne_pos,9);
    size_t first_class_occ_end = msg.find("::");
    size_t first_class_occ_beg = msg.rfind(" ",first_class_occ_end)+1;
    size_t length = first_class_occ_end-first_class_occ_beg+2;
    std::string class_name_scope = msg.substr(first_class_occ_beg,length);
    size_t class_name_scope_pos = std::string::npos;
    while ((class_name_scope_pos  = msg.find(class_name_scope,first_class_occ_end)) != std::string::npos)
        msg.erase(class_name_scope_pos,length);
    return msg;
}

LogScopeManager::LogScopeManager(std::string scope)
    : _scope(scope)
{
    Logger::increase_level(1);
    if ((!Logger::is_muted_at(0)) and Logger::prints_scope_entrance()) {
        Logger::println(0,"Enters '"+very_pretty_function(this->scope())+"'");
    }
}

std::string LogScopeManager::scope() const {
    return _scope;
}

LogScopeManager::~LogScopeManager() {
    if ((!Logger::is_muted_at(0)) and Logger::prints_scope_exit()) {
        Logger::println(0,"Exits '"+very_pretty_function(this->scope())+"'");
    }
    Logger::decrease_level(1);
    Logger::release(this->scope());
}

LoggerData::LoggerData(unsigned int current_level)
    : _current_level(current_level)
{ }

unsigned int LoggerData::current_level() const {
    return _current_level;
}

void LoggerData::enqueue_println(unsigned int level_increase, std::string text) {
    const std::lock_guard<std::mutex> lock(_queue_mutex);
    _raw_messages.push(LogRawMessage(std::string(),_current_level+level_increase,text));
}

void LoggerData::enqueue_hold(std::string scope, std::string text) {
    const std::lock_guard<std::mutex> lock(_queue_mutex);
    _raw_messages.push(LogRawMessage(scope,_current_level,text));
}

void LoggerData::enqueue_release(std::string scope) {
    const std::lock_guard<std::mutex> lock(_queue_mutex);
    _raw_messages.push(LogRawMessage(scope,_current_level,std::string()));
}

LogRawMessage LoggerData::dequeue() {
    const std::lock_guard<std::mutex> lock(_queue_mutex);
    LogRawMessage result = _raw_messages.front();
    _raw_messages.pop();
    return result;
}

void LoggerData::increase_level(unsigned int i) {
    _current_level += i;
}

void LoggerData::decrease_level(unsigned int i) {
    _current_level -= i;
}

unsigned int LoggerData::queue_size() const {
    return _raw_messages.size();
}

ImmediateLoggerScheduler::ImmediateLoggerScheduler() : _current_level(1) { }

unsigned int ImmediateLoggerScheduler::current_level() const {
    return _current_level;
}

void ImmediateLoggerScheduler::increase_level(unsigned int i) {
    _current_level+=i;
}

void ImmediateLoggerScheduler::decrease_level(unsigned int i) {
    _current_level-=i;
}

void ImmediateLoggerScheduler::println(unsigned int level_increase, std::string text) {
    Logger::_println(LogRawMessage(std::string(),_current_level+level_increase,text));
}

void ImmediateLoggerScheduler::hold(std::string scope, std::string text) {
    Logger::_hold(LogRawMessage(scope,_current_level,text));
}

void ImmediateLoggerScheduler::release(std::string scope) {
    Logger::_release(LogRawMessage(scope,_current_level,std::string()));
}

ConcurrentLoggerScheduler::ConcurrentLoggerScheduler() {
    _data.insert({std::this_thread::get_id(),SharedPointer<LoggerData>(new LoggerData(1))});
    _dequeueing_thread = std::thread(&ConcurrentLoggerScheduler::_dequeue_msgs, this);
    _terminate = false;
}

ConcurrentLoggerScheduler::~ConcurrentLoggerScheduler() {
    _terminate = true;
    _dequeueing_thread.join();
}

void ConcurrentLoggerScheduler::_create_data_instance(std::thread::id const& id) {
    _data.insert({id,SharedPointer<LoggerData>(new LoggerData(current_level()))});
}

SharedPointer<LoggerData> ConcurrentLoggerScheduler::_local_data() const {
    return _data.find(std::this_thread::get_id())->second;
}

unsigned int ConcurrentLoggerScheduler::current_level() const {
    return _local_data()->current_level();
}

void ConcurrentLoggerScheduler::increase_level(unsigned int i) {
    _local_data()->increase_level(i);
}

void ConcurrentLoggerScheduler::decrease_level(unsigned int i) {
    _local_data()->decrease_level(i);
}

void ConcurrentLoggerScheduler::println(unsigned int level_increase, std::string text) {
    _local_data()->enqueue_println(level_increase,text);
}

void ConcurrentLoggerScheduler::hold(std::string scope, std::string text) {
    _local_data()->enqueue_hold(scope,text);
}

void ConcurrentLoggerScheduler::release(std::string scope) {
    _local_data()->enqueue_release(scope);
}

std::pair<std::thread::id,unsigned int> ConcurrentLoggerScheduler::_largest_queue() {

    std::map<std::thread::id,SharedPointer<LoggerData>>::const_iterator it = _data.begin();

    std::pair<std::thread::id,unsigned int> result = std::make_pair(it->first,it->second->queue_size());
    for (++it; it != _data.end(); ++it) {
        if (it->second->queue_size() > result.second)
            result = std::make_pair(it->first,it->second->queue_size());
    }
    return result;
}

void ConcurrentLoggerScheduler::_dequeue_msgs() {
    while(true) {
        auto largest_queue = _largest_queue();
        if (largest_queue.second>0) {
            LogRawMessage msg = _data.find(largest_queue.first)->second->dequeue();

            switch (msg.kind()) {
                case RawMessageKind::PRINTLN :
                    Logger::_println(msg);
                    break;
                case RawMessageKind::HOLD :
                    Logger::_hold(msg);
                    break;
                case RawMessageKind::RELEASE :
                    Logger::_release(msg);
                    break;
                default:
                ARIADNE_FAIL_MSG("Unhandled RawMessageKind for printing.");
            }
        } else {
            if (not _terminate) std::this_thread::sleep_for(std::chrono::microseconds(200));
            else break;
        }
    }
}

void Logger::use_immediate_scheduler() {
    _scheduler.reset(new ImmediateLoggerScheduler());
}

void Logger::use_concurrent_scheduler() {
    _scheduler.reset(new ConcurrentLoggerScheduler());
}

void Logger::increase_level(unsigned int i) {
    _scheduler->increase_level(i);
}

void Logger::decrease_level(unsigned int i) {
    _scheduler->decrease_level(i);
}

void Logger::mute_increase_level() {
    _scheduler->increase_level(_MUTE_LEVEL_OFFSET);
}

void Logger::mute_decrease_level() {
    _scheduler->decrease_level(_MUTE_LEVEL_OFFSET);
}

bool Logger::is_muted_at(unsigned int i) {
    return (_verbosity < current_level()+i);
}

void Logger::set_verbosity(unsigned int v) {
    _verbosity = v;
}

void Logger::set_indents_based_on_level(bool b) {
    _indents_based_on_level = b;
}

void Logger::set_prints_level_on_change_only(bool b) {
    _prints_level_on_change_only = b;
}

void Logger::set_prints_scope_entrance(bool b) {
    _prints_scope_entrance = b;
}

void Logger::set_prints_scope_exit(bool b) {
    _prints_scope_exit = b;
}

void Logger::set_handles_multiline_output(bool b) {
    _handles_multiline_output = b;
}

void Logger::set_discards_newlines_and_indentation(bool b) {
    _discards_newlines_and_indentation = b;
}

unsigned int Logger::verbosity() {
    return _verbosity;
}

unsigned int Logger::current_level() {
    return _scheduler->current_level();
}

bool Logger::indents_based_on_level() {
    return _indents_based_on_level;
}

bool Logger::prints_level_on_change_only() {
    return _prints_level_on_change_only;
}

bool Logger::prints_scope_entrance() {
    return _prints_scope_entrance;
}

bool Logger::prints_scope_exit() {
    return _prints_scope_exit;
}

bool Logger::handles_multiline_output() {
    return _handles_multiline_output;
}

bool Logger::discards_newlines_and_indentation() {
    return _discards_newlines_and_indentation;
}

void Logger::println(unsigned int level_increase, std::string text) {
    _scheduler->println(level_increase, text);
}

void Logger::hold(std::string scope, std::string text) {
    _scheduler->hold(scope, text);
}

void Logger::release(std::string scope) {
    _scheduler->release(scope);
}

bool Logger::_is_holding() {
    return !_current_held_stack.empty();
}

unsigned int Logger::_get_window_columns() {
    struct winsize ws;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws);
    return ws.ws_col;
}

void Logger::_print_preamble_for_firstline(unsigned int level) {
    if (prints_level_on_change_only() and _cached_last_printed_level == level) std::clog << " ";
    else std::clog << level;
    std::clog << "| ";
    if (indents_based_on_level()) std::clog << std::string(level - 1, ' ');
}

void Logger::_print_preamble_for_extralines(unsigned int level) {
    std::clog << " · ";
    if (indents_based_on_level()) std::clog << std::string(level - 1, ' ');
}

std::string Logger::_discard_newlines_and_indentation(std::string const& text) {
    std::stringstream result;
    size_t text_ptr = 0;
    while(true) {
        std::size_t newline_pos = text.find('\n',text_ptr);
        if (newline_pos != std::string::npos) {
            size_t i=0;
            while(true) {
                if (text.at(newline_pos+1+i) == ' ') ++i;
                else break;
            }
            result << text.substr(text_ptr,newline_pos-text_ptr);
            text_ptr = newline_pos+1+i;
        } else {
            result << text.substr(text_ptr);
            break;
        }
    }
    return result.str();
}

void Logger::_print_held_line() {
    const unsigned int max_columns = _get_window_columns();
    unsigned int held_columns = 0;

    std::clog << '\r';
    for (auto msg : _current_held_stack) {
        held_columns = held_columns+(msg.level>9 ? 2 : 1)+3+msg.text.size();
        if (held_columns>max_columns+1) {
            std::string original = std::to_string(msg.level) + "| " + msg.text + " ";
            std::clog << original.substr(0,original.size()-(held_columns-max_columns+2)) << "..";
            held_columns=max_columns;
            break;
        } else if(held_columns==max_columns || held_columns==max_columns+1) {
            std::clog << msg.level << "| " << msg.text;
            held_columns=max_columns;
            break;
        } else {
            std::clog << msg.level << "| " << msg.text << " ";
        }
    }
    std::clog << std::flush;
    _cached_num_held_columns=held_columns;
}

void Logger::_println(LogRawMessage const& msg) {
    const unsigned int preamble_columns = 2+msg.level;
    if (_is_holding()) std::clog << '\r';

    _print_preamble_for_firstline(msg.level);

    std::string text = msg.text;
    if (discards_newlines_and_indentation()) text = _discard_newlines_and_indentation(text);

    if (handles_multiline_output() and msg.text.size() > 0) {
        const unsigned int max_columns = _get_window_columns();

        size_t text_ptr = 0;
        const size_t text_size = text.size();
        while(true) {
            if (text_size-text_ptr + preamble_columns > max_columns) {
                std::string to_print = text.substr(text_ptr,max_columns-preamble_columns);
                std::size_t newline_pos = to_print.find('\n');
                if (newline_pos != std::string::npos) {
                    std::clog << to_print.substr(0,newline_pos+1);
                    text_ptr += newline_pos+1;
                } else {
                    std::clog << to_print << '\n';
                    text_ptr += max_columns-preamble_columns;
                }
                _print_preamble_for_extralines(msg.level);
            } else {
                std::string to_print = text.substr(text_ptr,text_size-text_ptr);
                std::size_t newline_pos = to_print.find('\n');
                if (newline_pos != std::string::npos) {
                    std::clog << to_print.substr(0,newline_pos+1);
                    text_ptr += newline_pos+1;
                    _print_preamble_for_extralines(msg.level);
                } else {
                    std::clog << to_print;
                    break;
                }
            }
        }
    } else {
        std::clog << text;
    }

    if (_is_holding()) {
        const unsigned int printed_columns = preamble_columns+text.size();
        if (_cached_num_held_columns > printed_columns)
            std::clog << std::string(_cached_num_held_columns - printed_columns, ' ');
    }

    std::clog << '\n';
    if (_is_holding()) _print_held_line();

    _cached_last_printed_level = msg.level;
}

void Logger::_hold(LogRawMessage const& msg) {
    bool scope_found = false;
    for (unsigned int idx=0; idx<_current_held_stack.size(); ++idx) {
        if (_current_held_stack[idx].scope == msg.scope) { _current_held_stack[idx] = msg; scope_found = true; break; } }
    if (not scope_found) { _current_held_stack.push_back(msg); }
    _print_held_line();
}

void Logger::_release(LogRawMessage const& msg) {
    if (_is_holding()) {
        bool found = false;
        unsigned int i=0;
        for (; i<_current_held_stack.size(); ++i) {
            if (_current_held_stack[i].scope == msg.scope) {
                found = true;
                break;
            }
        }
        if (found) {
            unsigned int released_text_length = _current_held_stack[i].text.size()+4;
            std::vector<LogRawMessage> new_held_stack;
            if (i>0) {
                for (unsigned int j=0; j<i; ++j) {
                    new_held_stack.push_back(_current_held_stack[j]);
                }
            }
            if (i<_current_held_stack.size()-1) {
                for (unsigned int j=i+1; j<_current_held_stack.size(); ++j) {
                    new_held_stack.push_back(_current_held_stack[j]);
                }
            }
            _current_held_stack = new_held_stack;
            _print_held_line(); // Re-print
            std::clog << std::string(std::min(released_text_length,_get_window_columns()), ' '); // Fill the released chars with blanks
            if (not _is_holding()) // If nothing is held anymore, allow overwriting of the line
                std::clog << '\r';
            std::clog << std::flush;
        }
    }
}

} // namespace Ariadne

inline bool startup_function() {
    std::cerr << std::boolalpha;
    std::cout << std::boolalpha;
    std::clog << std::boolalpha;
    return true;
}

static const bool startup = startup_function();

