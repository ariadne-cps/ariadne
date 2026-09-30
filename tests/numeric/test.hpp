/***************************************************************************
 *            tests/numeric/test.hpp
 *
 *  Numeric test compatibility layer.
 *
 ****************************************************************************/

#ifndef ARIADNE_NUMERIC_TEST_HPP
#define ARIADNE_NUMERIC_TEST_HPP

#include "foundation/logical.hpp"
#include "utility/test.hpp"

#undef ARIADNE_TEST_ASSERT
#define ARIADNE_TEST_ASSERT(expression)                                 \
    {                                                                   \
        std::cout << #expression << ": " << std::flush;                 \
        auto result = (expression);                                     \
        if(definitely(result)) {                                        \
            std::cout << "true\n" << std::endl;                         \
        } else if(possibly(result)) {                                   \
            std::cout << "\nWARNING: indeterminate" << std::endl;       \
            std::cerr << "WARNING: " << __FILE__ << ":" << __LINE__ << ": " << ARIADNE_PRETTY_FUNCTION << ": Assertion `" << #expression << "' is indeterminate." << std::endl; \
        } else {                                                        \
            ++ARIADNE_TEST_FAILURES;                                    \
            std::cout << "\nERROR: false" << std::endl;                 \
            std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << ": " << ARIADNE_PRETTY_FUNCTION << ": Assertion `" << #expression << "' failed." << std::endl; \
        }                                                               \
    }

#undef ARIADNE_TEST_EQUALS
#define ARIADNE_TEST_EQUALS(expression,expected)                        \
    {                                                                   \
        std::cout << #expression << " == " << #expected << ": " << std::flush; \
        bool ok = decide((expression) == (expected));                   \
        if(ok) {                                                        \
            std::cout << "true\n" << std::endl;                         \
        } else {                                                        \
            ++ARIADNE_TEST_FAILURES;                                    \
            std::cout << "\nERROR: " << #expression << ":\n           " << (expression) << std::endl; \
            std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << ": " << ARIADNE_PRETTY_FUNCTION << ": Equality `" << #expression << " == " << #expected << "' failed;" << std::endl; \
            std::cerr << "  " << #expression << "=" << (expression) << std::endl; \
            std::cerr << "  " << #expected << "=" << (expected) << std::endl; \
        }                                                               \
    }

#undef ARIADNE_TEST_WITHIN
#define ARIADNE_TEST_WITHIN(expression,expected,tolerance)              \
    {                                                                   \
        std::cout << #expression << " ~ " << #expected << ": " << std::flush; \
        auto error = mag((expression) - (expected));                    \
        bool ok = decide(error <= (tolerance));                         \
        if(ok) {                                                        \
            std::cout << "true\n" << std::endl;                         \
        } else {                                                        \
            ++ARIADNE_TEST_FAILURES;                                    \
            std::cout << "\nERROR: " << #expression << ":\n           " << (expression) \
                      << "\n     : " << #expected << ":\n           " << (expected) \
                      << "\n     : error: " << error \
                      << "\n     : tolerance " << (tolerance) << std::endl; \
            std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << ": " << ARIADNE_PRETTY_FUNCTION << ": Approximate equality `" << #expression << " ~ " << #expected << "' failed." << std::endl; \
        }                                                               \
    }

#undef ARIADNE_TEST_BINARY_PREDICATE
#define ARIADNE_TEST_BINARY_PREDICATE(predicate,argument1,argument2)     \
    {                                                                   \
        std::cout << #predicate << "(" << (#argument1) << "," << (#argument2) << ") with " \
                  << #argument1 << "=" << (argument1) << ", " << #argument2 << "=" << (argument2) << ": " << std::flush; \
        bool ok = decide(predicate((argument1),(argument2)));           \
        if(ok) {                                                        \
            std::cout << "true\n" << std::endl;                         \
        } else {                                                        \
            ++ARIADNE_TEST_FAILURES;                                    \
            std::cout << "\nERROR: false" << std::endl;                 \
            std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << ": " << ARIADNE_PRETTY_FUNCTION << ": Predicate `" << #predicate << "(" << #argument1 << "," << #argument2 << ")' is false." << std::endl; \
        }                                                               \
    }


#undef ARIADNE_TEST_EQUAL
#define ARIADNE_TEST_EQUAL(expression1,expression2)                     \
    {                                                                   \
        std::cout << #expression1 << " == " << #expression2 << ": " << std::flush; \
        bool ok = decide((expression1) == (expression2));               \
        if(ok) {                                                        \
            std::cout << "true\n" << std::endl;                         \
        } else {                                                        \
            ++ARIADNE_TEST_FAILURES;                                    \
            std::cout << "\nERROR: " << #expression1 << ":\n           " << (expression1) \
                      << "\n     : " << #expression2 << ":\n           " << (expression2) << std::endl; \
            std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << ": " << ARIADNE_PRETTY_FUNCTION \
                      << ": Equality `" << #expression1 << " == " << #expression2 << "' failed." << std::endl; \
        }                                                               \
    }

#undef ARIADNE_TEST_UNARY_PREDICATE
#define ARIADNE_TEST_UNARY_PREDICATE(predicate,argument)                \
    {                                                                   \
        std::cout << #predicate << "(" << #argument << ") with " << #argument << "=" << (argument) << ": " << std::flush; \
        bool ok = decide(predicate((argument)));                        \
        if(ok) {                                                        \
            std::cout << "true\n" << std::endl;                         \
        } else {                                                        \
            ++ARIADNE_TEST_FAILURES;                                    \
            std::cout << "\nERROR: false" << std::endl;                 \
            std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << ": " << ARIADNE_PRETTY_FUNCTION \
                      << ": Predicate `" << #predicate << "(" << #argument << ")' is false." << std::endl; \
        }                                                               \
    }

#undef ARIADNE_TEST_COMPARE
#define ARIADNE_TEST_COMPARE(expression,comparison,expected)            \
    {                                                                   \
        std::cout << #expression << ": " << (expression) << std::flush; \
        bool ok = decide((expression) comparison (expected));           \
        if(ok) {                                                        \
            std::cout << " " << #comparison << " " << (expected) << ": true\n" << std::endl; \
        } else {                                                        \
            ++ARIADNE_TEST_FAILURES;                                    \
            std::cout << "\nERROR: expected: " << #expression << #comparison << #expected \
                      << "=" << (expected) << std::endl;                 \
            std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << ": " << ARIADNE_PRETTY_FUNCTION \
                      << ": Comparison `" << #expression << #comparison << #expected << "' failed; " \
                      << #expression << "=" << (expression) << "; " << #expected << "=" << (expected) << std::endl; \
        }                                                               \
    }

#endif
