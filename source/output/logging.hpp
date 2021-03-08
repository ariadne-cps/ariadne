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
#include "utility/macros.hpp"
#include "utility/pointer.hpp"
#include "utility/writable.hpp"

// Automatic level increase/decrease in a scope; meant to be used once within a function, at top scope; necessary for print holding.
#define ARIADNE_LOG_SCOPE_CREATE auto logscopemanager = LogScopeManager(ARIADNE_PRETTY_FUNCTION);
// Managed level increase/decrease around the function fn; if the function throws, manual decrease of the proper level is required.
#define ARIADNE_LOG_RUN_AT(level,fn) Logger::increase_level(level); fn; Logger::decrease_level(level);
// Mute the logger for the function fn; if the function throws, manual decrease of the proper level is required.
#define ARIADNE_LOG_RUN_MUTED(fn) Logger::mute_increase_level(); fn; Logger::mute_decrease_level();
// Print one line at the current level; the text shouldn't have carriage returns, but for efficiency purposes this is not checked.
#define ARIADNE_LOG_PRINTLN(text) { if (!Logger::is_muted_at(0)) { std::stringstream logger_stream; logger_stream << std::boolalpha << text; Logger::println(0,logger_stream.str()); } }
// Print one line at an increased level with respect to the current one; the text shouldn't have carriage returns, but for efficiency purposes this is not checked.
#define ARIADNE_LOG_PRINTLN_AT(level,text) { if (!Logger::is_muted_at(level)) { std::stringstream logger_stream; logger_stream << std::boolalpha << text; Logger::println(level,logger_stream.str()); } }
// Print a text at the bottom line, holding it until the function scope ends; this requires creation of the scope.
// Nested calls in separate functions append to the held line.
// The text for obvious reasons shouldn't have newlines and carriage returns; for efficiency purposes this is not checked.
#define ARIADNE_LOG_SCOPE_PRINTHOLD(text) { if (!Logger::is_muted_at(0)) { std::stringstream logger_stream; logger_stream << std::boolalpha << text; Logger::hold(logscopemanager.scope(),logger_stream.str()); } }

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

//! \brief A styling for a text character
//! \details Refer to https://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html#256-colors
//! for the detailed 256 map of colors for text and background
struct TerminalTextStyle {
    TerminalTextStyle(uint8_t fontcolor, uint8_t bgcolor, bool bold, bool underline);

    //! \brief Return the string code that applies the style
    std::string operator()() const;

    //! \brief If a style is actually set
    bool is_styled() const;

    uint8_t fontcolor; // The color of the character
    uint8_t bgcolor; // The color of the background of the character box
    bool bold; // If the character must be highlighted in bold
    bool underline; // If the character must be underlined

    inline static std::string RESET = "\u001b[0m";
};

static TerminalTextStyle TT_STYLE_NONE(0,0,false,false);
static TerminalTextStyle TT_STYLE_DARKBROWN(95,0,false,false);
static TerminalTextStyle TT_STYLE_DARKORANGE(130,0,false,false);
static TerminalTextStyle TT_STYLE_LIGHTBROWN(136,0,false,false);
static TerminalTextStyle TT_STYLE_BRIGHTORANGE(208,0,false,false);
static TerminalTextStyle TT_STYLE_CREAM(223,0,false,false);
static TerminalTextStyle TT_STYLE_OBSIDIAN(237,0,false,false);
static TerminalTextStyle TT_STYLE_DARKGREY(244,0,false,false);
static TerminalTextStyle TT_STYLE_LIGHTGREY(251,0,false,false);

//! \brief A set of text styles for different elements of the logging interface
struct TerminalTextTheme : public WritableInterface {
    TerminalTextTheme();
    TerminalTextTheme(TerminalTextStyle level_number, TerminalTextStyle level_shown_separator, TerminalTextStyle level_hidden_separator,
                      TerminalTextStyle multiline_separator, TerminalTextStyle assignment_comparison, TerminalTextStyle miscellaneous_operator,
                      TerminalTextStyle round_parentheses, TerminalTextStyle square_parentheses, TerminalTextStyle curly_parentheses,
                      TerminalTextStyle colon, TerminalTextStyle comma, TerminalTextStyle number, TerminalTextStyle keyword);

    TerminalTextStyle level_number; // The number representing the log level
    TerminalTextStyle level_shown_separator; // The separator when the level is shown
    TerminalTextStyle level_hidden_separator; // The separator when the level is hidden

    TerminalTextStyle multiline_separator; // The separator for multiple lines output
    TerminalTextStyle assignment_comparison; // Characters: =, >, <, !
    TerminalTextStyle miscellaneous_operator; // Characters: +, -, *, /, \, ^, |, &, %

    TerminalTextStyle round_parentheses; // Characters: ( )
    TerminalTextStyle square_parentheses; // Characters: [ ]
    TerminalTextStyle curly_parentheses; // Characters: { }

    TerminalTextStyle colon; // Character: :
    TerminalTextStyle comma; // Character: ,
    TerminalTextStyle number; // A 0-9 digit, or a dot after a digit, but not if parsing a numbered_variable

    TerminalTextStyle keyword; // Specific words, not pre/post-fixed by an alphanumeric character: virtual, const, true, false, inf

    bool has_style() const; // Whether at least one style is set

    virtual OutputStream& _write(OutputStream& os) const override;
};

//! \brief Empty theme, for not forcing any theme
static TerminalTextTheme TT_THEME_NONE = TerminalTextTheme();
//! \brief Theme for black background
static TerminalTextTheme TT_THEME_DARK = TerminalTextTheme(TT_STYLE_BRIGHTORANGE,TT_STYLE_BRIGHTORANGE,TT_STYLE_DARKGREY,
                                                           TT_STYLE_DARKGREY,TT_STYLE_BRIGHTORANGE,TT_STYLE_LIGHTBROWN,
                                                           TT_STYLE_DARKORANGE,TT_STYLE_DARKGREY,TT_STYLE_LIGHTBROWN,
                                                           TT_STYLE_LIGHTBROWN,TT_STYLE_DARKGREY,TT_STYLE_CREAM,
                                                           TT_STYLE_DARKGREY);
//! \brief Theme for white background
static TerminalTextTheme TT_THEME_LIGHT = TerminalTextTheme(TT_STYLE_BRIGHTORANGE,TT_STYLE_BRIGHTORANGE,TT_STYLE_DARKGREY,
                                                            TT_STYLE_DARKGREY,TT_STYLE_DARKBROWN,TT_STYLE_DARKBROWN,
                                                            TT_STYLE_DARKGREY,TT_STYLE_OBSIDIAN,TT_STYLE_LIGHTBROWN,
                                                            TT_STYLE_LIGHTBROWN,TT_STYLE_OBSIDIAN,TT_STYLE_DARKORANGE,
                                                            TT_STYLE_DARKGREY);

//! \brief Support class for managing log level increase/decrease in a scope
//! \details Since it is possible to capture the scope string of a function only, nested scopes in a function have the same scope string
//! (but work as expected in terms of level management)
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

    LogRawMessage(std::string scope, unsigned int level, std::string text);

    std::string scope;
    unsigned int level;
    std::string text;

    RawMessageKind kind() const;
};

class LoggerData;

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

//! \brief Configuration of visualisation settings for a Logger
class LoggerConfiguration : public WritableInterface {
  public:

    LoggerConfiguration();

    //! \brief Configuration setters

    //! \brief The depth of printing, where v=0 prevents all printing
    void set_verbosity(unsigned int v);
    //! \brief If true, indents each line with a number of spaces equal to the level
    void set_indents_based_on_level(bool b);
    //! \brief If false, the level is shown on each first line of a print; if true, only when a print has a different level
    void set_prints_level_on_change_only(bool b);
    //! \brief If true, each entrance in a defined scope is logged
    void set_prints_scope_entrance(bool b);
    //! \brief If true, each exit from a defined scope is logged
    void set_prints_scope_exit(bool b);
    //! \brief If true, text that would be printed on multiple lines is correctly preambled on each line with level and indentation
    void set_handles_multiline_output(bool b);
    //! \brief If true, discards from printing all newlines and their immediately following whitespaces
    void set_discards_newlines_and_indentation(bool b);

    //! \brief Configuration getters

    unsigned int verbosity() const;
    bool indents_based_on_level() const;
    bool prints_level_on_change_only() const;
    bool prints_scope_entrance() const;
    bool prints_scope_exit() const;
    bool handles_multiline_output() const;
    bool discards_newlines_and_indentation() const;

    //! \brief Style theme for terminal output
    void set_theme(TerminalTextTheme const& theme);
    //! \brief Get the current theme used
    TerminalTextTheme const& theme() const;
    //! \brief Add a keyword to the default ones offered in the them, forcing a given style
    //! \details Adding an existing keyword has no effect
    void add_custom_keyword(std::string const& text, TerminalTextStyle const& style);
    //! \brief Add a keyword to the default ones offered in the them, using the default style
    //! \details Adding an existing keyword has no effect
    void add_custom_keyword(std::string const& text);
    //! \brief Get the map of keywords (a TerminalTextStyle equal to TT_STYLE_NONE implies no custom style forced)
    std::map<std::string,TerminalTextStyle> const& custom_keywords() const;

    virtual OutputStream& _write(OutputStream& os) const override;

  private:
    unsigned int _verbosity;
    bool _indents_based_on_level;
    bool _prints_level_on_change_only;
    bool _prints_scope_entrance;
    bool _prints_scope_exit;
    bool _handles_multiline_output;
    bool _discards_newlines_and_indentation;
    TerminalTextTheme _theme;
    std::map<std::string,TerminalTextStyle> _custom_keywords;
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
    static unsigned int cached_last_printed_level();

    static LoggerConfiguration& configuration();

  private:
    static bool _is_assignment_comparison(char const& c);
    static std::string _apply_theme(std::string const& text);
    static std::string _apply_theme_for_keywords(std::string const& text);
    static void _print_preamble_for_firstline(unsigned int level);
    static void _print_preamble_for_extralines(unsigned int level);
    static std::string _discard_newlines_and_indentation(std::string const& text);
    static void _cover_held_columns_with_whitespaces(unsigned int printed_columns);
    static void _print_held_line();
    static void _println(LogRawMessage const& msg);
    static void _hold(LogRawMessage const& msg);
    static void _release(LogRawMessage const& msg);
    static unsigned int _get_window_columns();
    static bool _is_holding();
  private:
    const static unsigned int _MUTE_LEVEL_OFFSET = 1024;
    inline static std::vector<LogRawMessage> _current_held_stack;
    inline static unsigned int _cached_num_held_columns = 0;
    inline static unsigned int _cached_last_printed_level = 0;
    inline static SharedPointer<LoggerSchedulerInterface> _scheduler = SharedPointer<LoggerSchedulerInterface>(new ConcurrentLoggerScheduler());
    inline static LoggerConfiguration _configuration;
};

// Global log output file
extern std::ofstream log_file_stream;

//! \brief Redirect logging output to file \a filename.
void redirect_log(const char* filename);

} // namespace Ariadne

#endif // ARIADNE_LOGGING_HPP
