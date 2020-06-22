/***************************************************************************
 *            output/logging.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
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

/*! \file output/logging.hpp
 *  \brief Support for writing debugging output to a logging stream.
 */

#ifndef ARIADNE_LOGGING_HPP
#define ARIADNE_LOGGING_HPP

#include <iostream>
#include <fstream>
#include <cstring>
#include <utility>
#include <vector>
#include <queue>
#include <map>
#include <sstream>
#include <thread>
#include <mutex>
#include <atomic>
#include "utility/pointer.hpp"

// Automatic level increase/decrease in a scope; meant to be used once within a function, at top scope; necessary for print holding.
#define ARIADNE_LOG_SCOPE_CREATE auto logscopemanager = LogScopeManager(ARIADNE_PRETTY_FUNCTION);
// Managed level increase/decrease around the function fn; if the function throws, manual decrease of the proper level is required.
#define ARIADNE_LOG_RUN_AT(level,fn) Logger::increase_level(level); fn; Logger::decrease_level(level);
// Mute the logger for the function fn; if the function throws, manual decrease of the proper level is required.
#define ARIADNE_LOG_RUN_MUTED(fn) Logger::mute_increase_level(); fn; Logger::mute_decrease_level();
// Print one line at the current level; the text shouldn't have carriage returns, but for efficiency purposes this is not checked.
#define ARIADNE_LOG_PRINTLN(text) { if (!Logger::is_muted_at(0)) { std::stringstream logger_stream; logger_stream << text; Logger::println(0,logger_stream.str()); } }
// Print one line at an increased level with respect to the current one; the text shouldn't have carriage returns, but for efficiency purposes this is not checked.
#define ARIADNE_LOG_PRINTLN_AT(level,text) { if (!Logger::is_muted_at(level)) { std::stringstream logger_stream; logger_stream << text; Logger::println(level,logger_stream.str()); } }
// Print a text at the bottom line, holding it until the function scope ends; this requires creation of the scope.
// Nested calls in separate functions append to the held line.
// The text for obvious reasons shouldn't have newlines and carriage returns; for efficiency purposes this is not checked.
#define ARIADNE_LOG_SCOPE_PRINTHOLD(text) { if (!Logger::is_muted_at(0)) { std::stringstream logger_stream; logger_stream << text; Logger::hold(logscopemanager.scope(),logger_stream.str()); } }

namespace Ariadne {

inline unsigned int get_verbosity(int argc, const char* argv[]) {
    if(argc>1) {
        if(std::strcmp(argv[1],"-v")==0) {
            if(argc>2) {
                int val = std::atoi(argv[2]);
                if (val < 0) std::cerr << "Verbosity should be a non-negative value.\n";
                return static_cast<unsigned int>(val);
            }
        } else {
            std::cerr << "Unrecognised command-line option \"" << argv[1] << "\"\n";
        }
    }
    return 0;
}


class LogScopeManager {
  public:
    LogScopeManager(std::string scope);
    std::string scope() const;
  public:
    virtual ~LogScopeManager();
  private:
    mutable std::string _scope;
};

//! \brief Supported kinds of log messages
enum RawMessageKind { PRINTLN, HOLD, RELEASE };

//! \brief A log message in raw form, as submitted, before formatting for actual output
struct LogRawMessage {

    LogRawMessage(std::string scope_, unsigned int level_, std::string text_) : scope(scope_), level(level_), text(text_) { }

    std::string scope;
    unsigned int level;
    std::string text;

    RawMessageKind kind() {
        if (scope.empty()) return RawMessageKind::PRINTLN;
        else if (!text.empty()) return RawMessageKind::HOLD;
        else return RawMessageKind::RELEASE;
    }
};

class ConcurrentLoggerScheduler;

//! \brief Thread-based enqueued log data
class LoggerData {
    friend class ConcurrentLoggerScheduler;
protected:
    LoggerData(unsigned int current_level);

    void enqueue_println(unsigned int level_increase, std::string text);
    void enqueue_hold(std::string scope, std::string text);
    void enqueue_release(std::string scope);

    LogRawMessage dequeue();

    void increase_level(unsigned int i);
    void decrease_level(unsigned int i);

    unsigned int current_level() const;

    unsigned int queue_size() const;
private:
    unsigned int _current_level;
    std::queue<LogRawMessage> _raw_messages;
    std::mutex _queue_mutex;
};

class LoggerSchedulerInterface {
  public:
    virtual void println(unsigned int level_increase, std::string text) = 0;
    virtual void hold(std::string scope, std::string text) = 0;
    virtual void release(std::string scope) = 0;
    virtual unsigned int current_level() const = 0;
    virtual void increase_level(unsigned int i) = 0;
    virtual void decrease_level(unsigned int i) = 0;
    virtual ~LoggerSchedulerInterface() = default;
};

//! \brief A Logger scheduler that prints immediately. Not designed for concurrency, since
//! the current level should change based on the thread. Also clearly the outputs can overlap arbitrarily.
class ImmediateLoggerScheduler : public LoggerSchedulerInterface {
public:
    ImmediateLoggerScheduler();
    virtual void println(unsigned int level_increase, std::string text) override;
    virtual void hold(std::string scope, std::string text) override;
    virtual void release(std::string scope) override;
    virtual unsigned int current_level() const override;
    virtual void increase_level(unsigned int i) override;
    virtual void decrease_level(unsigned int i) override;
private:
    unsigned int _current_level;
};

//! \brief A Logger scheduler that enqueues messages and prints them sequentially.
//! The order of printing takes into account the queue for each thread.
class ConcurrentLoggerScheduler : public LoggerSchedulerInterface {
  public:
    ConcurrentLoggerScheduler();
    virtual void println(unsigned int level_increase, std::string text) override;
    virtual void hold(std::string scope, std::string text) override;
    virtual void release(std::string scope) override;
    virtual unsigned int current_level() const override;
    virtual void increase_level(unsigned int i) override;
    virtual void decrease_level(unsigned int i) override;
    ~ConcurrentLoggerScheduler();
  private:
    SharedPointer<LoggerData> _local_data() const;
    std::pair<std::thread::id,unsigned int> _largest_queue();
    void _dequeue_msgs();
    void _create_data_instance(std::thread::id const& id);
  private:
    std::map<std::thread::id,SharedPointer<LoggerData>> _data;
    std::thread _dequeueing_thread;
    std::atomic<bool> _terminate;
};

//! \brief A static class for log output handling.
//! Configuration and final printing is done here, while scheduling is
//! done in a separate class.
class Logger {
    friend class ConcurrentLoggerScheduler;
    friend class ImmediateLoggerScheduler;
  public:
    static void use_immediate_scheduler();
    static void use_concurrent_scheduler();

    static void println(unsigned int level_increase, std::string text);
    static void hold(std::string scope, std::string text);
    static void release(std::string scope);

    static void increase_level(unsigned int i);
    static void decrease_level(unsigned int i);
    static void mute_increase_level();
    static void mute_decrease_level();

    static bool is_muted_at(unsigned int i);

    static unsigned int current_level();

  public:
    //! \brief Configuration setters

    //! \brief The depth of printing, where v=0 prevents all printing
    static void set_verbosity(unsigned int v);
    //! \brief If true, indents each line with a number of spaces equal to the level
    static void set_indents_based_on_level(bool b);
    //! \brief If false, the level is shown on each first line of a print; if true, only when a print has a different level
    static void set_prints_level_on_change_only(bool b);
    //! \brief If true, each entrance in a defined scope is logged
    static void set_prints_scope_entrance(bool b);
    //! \brief If true, each exit from a defined scope is logged
    static void set_prints_scope_exit(bool b);
    //! \brief If true, text that would be printed on multiple lines is correctly preambled on each line with level and indentation
    static void set_handles_multiline_output(bool b);
    //! \brief If true, discards from printing all newlines and their immediately following whitespaces
    static void set_discards_newlines_and_indentation(bool b);
    // Configuration getters
    static unsigned int verbosity();
    static bool indents_based_on_level();
    static bool prints_level_on_change_only();
    static bool prints_scope_entrance();
    static bool prints_scope_exit();
    static bool handles_multiline_output();
    static bool discards_newlines_and_indentation();
  private:
    static void _print_preamble_for_firstline(unsigned int level);
    static void _print_preamble_for_extralines(unsigned int level);
    static std::string _discard_newlines_and_indentation(std::string const& text);
    static void _print_held_line();
    static void _println(LogRawMessage const& msg);
    static void _hold(LogRawMessage const& msg);
    static void _release(LogRawMessage const& msg);
    static unsigned int _get_window_columns();
    static bool _is_holding();

  private:
    const static unsigned int _MUTE_LEVEL_OFFSET = 1024;

    inline static unsigned int _verbosity = 0;
    inline static bool _indents_based_on_level = true;
    inline static bool _prints_level_on_change_only = true;
    inline static bool _prints_scope_entrance = false;
    inline static bool _prints_scope_exit = false;
    inline static bool _handles_multiline_output = true;
    inline static bool _discards_newlines_and_indentation = false;
    inline static std::vector<LogRawMessage> _current_held_stack;
    inline static unsigned int _cached_num_held_columns = 0;
    inline static unsigned int _cached_last_printed_level = 0;
    inline static SharedPointer<LoggerSchedulerInterface> _scheduler = SharedPointer<LoggerSchedulerInterface>(new ConcurrentLoggerScheduler());
};

// Global log output file
extern std::ofstream log_file_stream;

//! \brief Redirect logging output to file \a filename.
void redirect_log(const char* filename);

} // namespace Ariadne

#endif // ARIADNE_LOGGING_HPP
