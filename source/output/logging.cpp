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

#include "config.hpp"
#include "output/logging.hpp"
#include "utility/writable.hpp"
#include "utility/macros.hpp"

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

TerminalTextStyle::TerminalTextStyle(uint8_t fontcolor_, uint8_t bgcolor_, bool bold_, bool underline_) :
    fontcolor(fontcolor_), bgcolor(bgcolor_), bold(bold_), underline(underline_)
{ }

std::string TerminalTextStyle::operator()() const {
    std::stringstream ss;
    if (((int)fontcolor) > 0) ss << "\u001b[38;5;" << (int)fontcolor << "m";
    if (((int)bgcolor) > 0) ss << "\u001b[48;5;" << (int)fontcolor << "m";
    if (bold) ss << "\u001b[1m";
    if (underline) ss << "\u001b[4m";
    return ss.str();
}

bool TerminalTextStyle::is_styled() const {
    return ((int)fontcolor) > 0 or ((int)bgcolor) > 0 or bold or underline;
}

TerminalTextTheme::TerminalTextTheme(TerminalTextStyle level_number_, TerminalTextStyle level_shown_separator_, TerminalTextStyle level_hidden_separator_,
                                     TerminalTextStyle multiline_separator_, TerminalTextStyle assignment_comparison_, TerminalTextStyle miscellaneous_operator_,
                                     TerminalTextStyle round_parentheses_, TerminalTextStyle square_parentheses_, TerminalTextStyle curly_parentheses_,
                                     TerminalTextStyle colon_, TerminalTextStyle comma_, TerminalTextStyle number_, TerminalTextStyle keyword_) :
        level_number(level_number_), level_shown_separator(level_shown_separator_), level_hidden_separator(level_hidden_separator_),
        multiline_separator(multiline_separator_), assignment_comparison(assignment_comparison_), miscellaneous_operator(miscellaneous_operator_),
        round_parentheses(round_parentheses_), square_parentheses(square_parentheses_), curly_parentheses(curly_parentheses_),
        colon(colon_), comma(comma_), number(number_), keyword(keyword_) { }

TerminalTextTheme::TerminalTextTheme() :
        level_number(TT_STYLE_NONE), level_shown_separator(TT_STYLE_NONE), level_hidden_separator(TT_STYLE_NONE),
        multiline_separator(TT_STYLE_NONE), assignment_comparison(TT_STYLE_NONE), miscellaneous_operator(TT_STYLE_NONE),
        round_parentheses(TT_STYLE_NONE), square_parentheses(TT_STYLE_NONE), curly_parentheses(TT_STYLE_NONE),
        colon(TT_STYLE_NONE), comma(TT_STYLE_NONE), number(TT_STYLE_NONE), keyword(TT_STYLE_NONE) { }

bool TerminalTextTheme::has_style() const {
    return (level_number.is_styled() or level_shown_separator.is_styled() or level_hidden_separator.is_styled() or
            multiline_separator.is_styled() or assignment_comparison.is_styled() or miscellaneous_operator.is_styled() or
            round_parentheses.is_styled() or square_parentheses.is_styled() or curly_parentheses.is_styled() or
            colon.is_styled() or comma.is_styled() or number.is_styled() or keyword.is_styled());
}

OutputStream& TerminalTextTheme::_write(OutputStream& os) const {
    // This should not be used by ARIADNE_LOG_PRINTLN, otherwise it will be parsed further, breaking level separator styling
    os << "TerminalTextTheme("
       << "\n  level_number= " << level_number() << "1 2 3 4 5 6 7 8 9" << TerminalTextStyle::RESET
            << ",\n  level_shown_separator= " << level_shown_separator() << "|" << TerminalTextStyle::RESET
            << ",\n  level_hidden_separator= " << level_hidden_separator() << "|" << TerminalTextStyle::RESET
            << ",\n  multiline_separator= " << multiline_separator() << "·" << TerminalTextStyle::RESET
            << ",\n  assignment_comparison= " << assignment_comparison() << "= > < !" << TerminalTextStyle::RESET
            << ",\n  miscellaneous_operator= " << miscellaneous_operator() << "+ - * / \\ ^ | & %" << TerminalTextStyle::RESET
            << ",\n  round_parentheses= " << round_parentheses() << "( )" << TerminalTextStyle::RESET
            << ",\n  square_parentheses= " << square_parentheses() << "[ ]" << TerminalTextStyle::RESET
            << ",\n  curly_parentheses= " << curly_parentheses() << "{ }" << TerminalTextStyle::RESET
            << ",\n  comma= " << comma() << ":" << TerminalTextStyle::RESET
            << ",\n  number= " << number() << "1.2" << TerminalTextStyle::RESET
            << ",\n  keyword= " << keyword() << "virtual const true false inf" << TerminalTextStyle::RESET
       << "\n)";
    return os;
}

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

LogScopeManager::LogScopeManager(std::string scope)
    : _scope(scope)
{
    Logger::increase_level(1);
    if ((!Logger::is_muted_at(0)) and Logger::configuration().prints_scope_entrance()) {
        Logger::println(0,"Enters '"+very_pretty_function(this->scope())+"'");
    }
}

std::string LogScopeManager::scope() const {
    return _scope;
}

LogScopeManager::~LogScopeManager() {
    if ((!Logger::is_muted_at(0)) and Logger::configuration().prints_scope_exit()) {
        Logger::println(0,"Exits '"+very_pretty_function(this->scope())+"'");
    }
    Logger::decrease_level(1);
    Logger::release(this->scope());
}

LogRawMessage::LogRawMessage(std::string scope_, unsigned int level_, std::string text_) :
    scope(scope_), level(level_), text(text_)
{ }

RawMessageKind LogRawMessage::kind() const {
    if (scope.empty()) return RawMessageKind::PRINTLN;
    else if (!text.empty()) return RawMessageKind::HOLD;
    else return RawMessageKind::RELEASE;
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
            // Allow the thread to sleep, for a larger time in the case of high verbosity (for verbosity 1, a 50ms sleep)
            // The current level is used as an average of the levels printed
            if (not _terminate) std::this_thread::sleep_for(std::chrono::microseconds(std::max(1,100000>>Logger::cached_last_printed_level())));
            else break;
        }
    }
}

LoggerConfiguration::LoggerConfiguration() :
    _verbosity(0), _indents_based_on_level(true), _prints_level_on_change_only(true), _prints_scope_entrance(false),
    _prints_scope_exit(false), _handles_multiline_output(true), _discards_newlines_and_indentation(false),
    _theme(TT_THEME_NONE)
{ }

LoggerConfiguration& Logger::configuration() {
    return _configuration;
}

void LoggerConfiguration::set_verbosity(unsigned int v) {
    _verbosity = v;
}

void LoggerConfiguration::set_indents_based_on_level(bool b) {
    _indents_based_on_level = b;
}

void LoggerConfiguration::set_prints_level_on_change_only(bool b) {
    _prints_level_on_change_only = b;
}

void LoggerConfiguration::set_prints_scope_entrance(bool b) {
    _prints_scope_entrance = b;
}

void LoggerConfiguration::set_prints_scope_exit(bool b) {
    _prints_scope_exit = b;
}

void LoggerConfiguration::set_handles_multiline_output(bool b) {
    _handles_multiline_output = b;
}

void LoggerConfiguration::set_discards_newlines_and_indentation(bool b) {
    _discards_newlines_and_indentation = b;
}

void LoggerConfiguration::set_theme(TerminalTextTheme const& theme) {
    _theme = theme;
}

unsigned int LoggerConfiguration::verbosity() const {
    return _verbosity;
}

bool LoggerConfiguration::indents_based_on_level() const {
    return _indents_based_on_level;
}

bool LoggerConfiguration::prints_level_on_change_only() const {
    return _prints_level_on_change_only;
}

bool LoggerConfiguration::prints_scope_entrance() const {
    return _prints_scope_entrance;
}

bool LoggerConfiguration::prints_scope_exit() const {
    return _prints_scope_exit;
}

bool LoggerConfiguration::handles_multiline_output() const {
    return _handles_multiline_output;
}

bool LoggerConfiguration::discards_newlines_and_indentation() const {
    return _discards_newlines_and_indentation;
}

TerminalTextTheme const& LoggerConfiguration::theme() const {
    return _theme;
}

void LoggerConfiguration::add_custom_keyword(std::string const& text, TerminalTextStyle const& style) {
    _custom_keywords.insert({text,style});
}

void LoggerConfiguration::add_custom_keyword(std::string const& text) {
    add_custom_keyword(text,TT_STYLE_NONE);
}

std::map<std::string,TerminalTextStyle> const& LoggerConfiguration::custom_keywords() const {
    return _custom_keywords;
}

OutputStream& LoggerConfiguration::_write(OutputStream& os) const {
    os << "LoggerConfiguration("
       << "\n  verbosity=" << _verbosity
       << ",\n  indents_based_on_level=" << _indents_based_on_level
       << ",\n  prints_level_on_change_only=" << _prints_level_on_change_only
       << ",\n  prints_scope_entrance=" << _prints_scope_entrance
       << ",\n  prints_scope_exit=" << _prints_scope_exit
       << ",\n  handles_multiline_output=" << _handles_multiline_output
       << ",\n  discards_newlines_and_indentation=" << _discards_newlines_and_indentation
       << ",\n  theme=(not shown)" // To show theme colors appropriately, print the theme object directly on standard output
       << "\n)";
    return os;
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
    return (_configuration.verbosity() < current_level()+i);
}

unsigned int Logger::current_level() {
    return _scheduler->current_level();
}

unsigned int Logger::cached_last_printed_level() {
    return _cached_last_printed_level;
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

std::string Logger::_apply_theme(std::string const& text) {
    auto theme = _configuration.theme();
    if (theme.has_style()) {
        std::stringstream ss;
        for(auto it = text.begin(); it != text.end(); ++it) {
            const char& c = *it;
            switch (c) {
                case '=':
                case '>':
                case '<':
                case '!':
                    ss << theme.assignment_comparison() << c << TerminalTextStyle::RESET;
                    break;
                case '(':
                case ')':
                    ss << theme.round_parentheses() << c << TerminalTextStyle::RESET;
                    break;
                case '[':
                case ']':
                    ss << theme.square_parentheses() << c << TerminalTextStyle::RESET;
                    break;
                case '{':
                case '}':
                    ss << theme.curly_parentheses() << c << TerminalTextStyle::RESET;
                    break;
                case ':':
                    ss << theme.colon() << c << TerminalTextStyle::RESET;
                    break;
                case ',':
                    ss << theme.comma() << c << TerminalTextStyle::RESET;
                    break;
                case '.':
                    if (it != text.begin() and isdigit(*(it-1)))
                        ss << theme.number() << c << TerminalTextStyle::RESET;
                    else
                        ss << c;
                    break;
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                {
                    // Exclude strings that end with a number (supported up to 2 digits) to account for numbered variables
                    // For simplicity, this does not work across multiple lines
                    bool styled = true;
                    if (it != text.begin()) {
                        if (isalpha(*(it - 1)))
                            styled = false;
                        else if (isdigit(*(it - 1))) {
                            if ((it - 1) != text.begin()) {
                                if (isalpha(*(it - 2)))
                                    styled = false;
                            }
                        }
                    }
                    if (styled) ss << theme.number() << c << TerminalTextStyle::RESET;
                    else ss << c;
                    break;
                }
                case '+':
                case '-':
                case '*':
                case '/':
                case '\\':
                case '^':
                case '|':
                case '&':
                case '%':
                    ss << theme.miscellaneous_operator() << c << TerminalTextStyle::RESET;
                    break;
                default:
                    ss << c;
            }
        }
        return _apply_theme_for_keywords(ss.str());
    } else return text;
}

bool isalpha_withstylecodes(std::string text, size_t pos) {
    auto c = text.at(pos);
    if (not isalpha(c)) {
        return false;
    } else if (c == 'm' and pos > 2) { // If this is the last character of a feasible style code
        auto sub = text.substr(pos-3,3);
        // Check for reset or bold or underline
        if (sub == "\u001b[0" or sub == "\u001b[1" or sub == "\u001b[4") return false;
        else if (isdigit(text.at(pos-1)) and pos > 8) { // If this is a feasible font/bg color (at 1,2 or 3 digits)
            if (text.substr(pos-8,7) == "\u001b[38;5;" or text.substr(pos-8,7) == "\u001b[48;5;") return false;
            else if (isdigit(text.at(pos-2)) and pos > 9) {
                if (text.substr(pos-9,7) == "\u001b[38;5;" or text.substr(pos-9,7) == "\u001b[48;5;") return false;
                else if (isdigit(text.at(pos-3)) and pos > 10) {
                    if (text.substr(pos-10,7) == "\u001b[38;5;" or text.substr(pos-10,7) == "\u001b[48;5;") return false;
                    else return true;
                } else return true;
            } else return true;
        } else return true;
    } else return true;
}

std::string Logger::_apply_theme_for_keywords(std::string const& text) {
    std::string result = text;
    std::map<std::string,TerminalTextStyle> keyword_styles = {{"virtual",_configuration.theme().keyword},
                                                              {"const",_configuration.theme().keyword},
                                                              {"true",_configuration.theme().keyword},
                                                              {"false",_configuration.theme().keyword},
                                                              {"inf",_configuration.theme().keyword}};
    keyword_styles.insert(_configuration.custom_keywords().begin(),_configuration.custom_keywords().end());

    for (auto kws : keyword_styles) {
        std::stringstream current_result;
        size_t kw_length = kws.first.length();
        size_t kw_pos = std::string::npos;
        size_t scan_pos = 0;
        while ((kw_pos  = result.find(kws.first,scan_pos)) != std::string::npos) {
            // Apply the theme only if the keyword is not adjacent to an alphanumerical character
            if ((kw_pos > 0 and (isalpha_withstylecodes(result,kw_pos-1) or isdigit(result.at(kw_pos-1)))) or
                (kw_pos+kw_length < result.length() and (isalpha(result.at(kw_pos+kw_length)) or isdigit(result.at(kw_pos+kw_length)))))
                current_result << result.substr(scan_pos,kw_pos-scan_pos+kw_length);
            else
                current_result << result.substr(scan_pos,kw_pos-scan_pos) << kws.second() << result.substr(kw_pos,kw_length) << TerminalTextStyle::RESET;

            scan_pos = kw_pos+kw_length;
        }
        if (scan_pos != 0) {
            if (scan_pos < result.length()) current_result << result.substr(scan_pos);
            result = current_result.str();
        }
    }
    return result;
}

void Logger::_print_preamble_for_firstline(unsigned int level) {
    auto theme = _configuration.theme();
    bool level_printed = true;
    if (_configuration.prints_level_on_change_only() and _cached_last_printed_level == level) {
        level_printed = false;
        std::clog << " ";
    }
    else {
        if (theme.level_number.is_styled()) std::clog << theme.level_number() << level << TerminalTextStyle::RESET;
        else std::clog << level;
    }
    if (not level_printed and theme.level_hidden_separator.is_styled()) {
        std::clog << theme.level_hidden_separator() << "|" << TerminalTextStyle::RESET;
    } else if (level_printed and theme.level_shown_separator.is_styled()) {
        std::clog << theme.level_shown_separator() << "|" << TerminalTextStyle::RESET;
    } else {
        std::clog << "|";
    }
    if (_configuration.indents_based_on_level()) std::clog << std::string(level, ' ');
}

void Logger::_print_preamble_for_extralines(unsigned int level) {
    auto theme = _configuration.theme();
    std::clog << " ";
    if (theme.multiline_separator.is_styled()) std::clog << theme.multiline_separator() << "·" << TerminalTextStyle::RESET;
    else std::clog << "·";

    if (_configuration.indents_based_on_level()) std::clog << std::string(level, ' ');
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
    auto theme = _configuration.theme();
    const unsigned int max_columns = _get_window_columns();
    unsigned int held_columns = 0;

    std::clog << '\r';
    for (auto msg : _current_held_stack) {
        held_columns = held_columns+(msg.level>9 ? 2 : 1)+3+msg.text.size();
        if (held_columns>max_columns+1) {
            std::string original = theme.level_number() + std::to_string(msg.level) + TerminalTextStyle::RESET +
                                   theme.level_shown_separator() + "|" + TerminalTextStyle::RESET + " " + _apply_theme(msg.text) + " ";
            std::clog << original.substr(0,original.size()-(held_columns-max_columns+2)) << "..";
            held_columns=max_columns;
            break;
        } else if(held_columns==max_columns || held_columns==max_columns+1) {
            std::clog << theme.level_number() << msg.level << TerminalTextStyle::RESET <<
                         theme.level_shown_separator() << "|" << TerminalTextStyle::RESET << " " << _apply_theme(msg.text);
            held_columns=max_columns;
            break;
        } else {
            std::clog << theme.level_number() << msg.level << TerminalTextStyle::RESET <<
                         theme.level_shown_separator() << "|" << TerminalTextStyle::RESET << " " << _apply_theme(msg.text) << " ";
        }
    }
    std::clog << std::flush;
    _cached_num_held_columns=held_columns;
    // Sleep for a time exponential in the last printed level (used as an average of the levels printed)
    // In this way, holding is printed more cleanly against the OS buffering, which overrules flushing
    std::this_thread::sleep_for(std::chrono::microseconds(10<<_cached_last_printed_level));
}

void Logger::_cover_held_columns_with_whitespaces(unsigned int printed_columns) {
    if (_is_holding()) {
        if (_cached_num_held_columns > printed_columns)
            std::clog << std::string(_cached_num_held_columns - printed_columns, ' ');
    }
}

void Logger::_println(LogRawMessage const& msg) {
    const unsigned int preamble_columns = (msg.level>9 ? 3:2)+msg.level;

    // If holding, we must write over the held line first
    if (_is_holding()) std::clog << '\r';

    _print_preamble_for_firstline(msg.level);

    std::string text = msg.text;
    if (configuration().discards_newlines_and_indentation()) text = _discard_newlines_and_indentation(text);

    if (configuration().handles_multiline_output() and msg.text.size() > 0) {
        const unsigned int max_columns = _get_window_columns();

        size_t text_ptr = 0;
        const size_t text_size = text.size();
        while(true) {
            if (text_size-text_ptr + preamble_columns > max_columns) { // (remaining) Text too long for a single terminal line
                std::string to_print = text.substr(text_ptr,max_columns-preamble_columns);
                std::size_t newline_pos = to_print.find('\n');
                if (newline_pos != std::string::npos) { // A newline is found before reaching the end of the terminal line
                    std::clog << _apply_theme(to_print.substr(0,newline_pos));
                    _cover_held_columns_with_whitespaces(preamble_columns+to_print.substr(0,newline_pos).size());
                    text_ptr += newline_pos+1;
                } else { // Text reaches the end of the terminal line
                    std::clog << _apply_theme(to_print);
                    _cover_held_columns_with_whitespaces(preamble_columns+to_print.size());
                    text_ptr += max_columns-preamble_columns;
                }
                std::clog << '\n';
                if (_is_holding()) _print_held_line();
                if (_is_holding()) std::clog << '\r';

                _print_preamble_for_extralines(msg.level);
            } else { // (remaining) Text shorter than the terminal line
                std::string to_print = text.substr(text_ptr,text_size-text_ptr);
                std::size_t newline_pos = to_print.find('\n');
                if (newline_pos != std::string::npos) { // A newline is found before reaching the end of the terminal line
                    std::clog << _apply_theme(to_print.substr(0,newline_pos));
                    _cover_held_columns_with_whitespaces(preamble_columns+to_print.substr(0,newline_pos).size());
                    std::clog << '\n';
                    if (_is_holding()) {
                        _print_held_line();
                        std::clog << '\r';
                    }

                    text_ptr += newline_pos+1;
                    _print_preamble_for_extralines(msg.level);
                } else { // Text reaches the end of the terminal line
                    std::clog << _apply_theme(to_print);
                    _cover_held_columns_with_whitespaces(preamble_columns+to_print.size());
                    std::clog << '\n';
                    if (_is_holding()) _print_held_line();

                    break;
                }
            }
        }
    } else { // No multiline is handled, \n characters are handled by the terminal
        std::clog << _apply_theme(text);
        _cover_held_columns_with_whitespaces(preamble_columns+text.size());
        std::clog << '\n';
        if (_is_holding()) _print_held_line();
    }

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

