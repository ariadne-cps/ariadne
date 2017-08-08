/***************************************************************************
 *            utility.cpp
 *
 *  Copyright 2017  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "stack_trace.hpp"

#include <cassert>
#include <iostream>
#include <execinfo.h>
#include <dlfcn.h>
#include <cxxabi.h>

namespace Ariadne {

void stack_trace() {
    static const unsigned int CALLSTACK_SIZE = 128;
    static const unsigned int BACKTRACE_BUFFER_SIZE = 1024;
    int skip=0;
    void *callstack[CALLSTACK_SIZE];
    const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
    char buffer[BACKTRACE_BUFFER_SIZE];
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);
    assert (symbols != nullptr);
//    std::cerr << "\nnFrames=" << nFrames << "\n";
    for (int i = skip; i < nFrames; i++) {
        Dl_info info;
        if (dladdr(callstack[i], &info)) {
            const char* dli_sname = info.dli_sname;
            char* no_output_buffer = nullptr;
            size_t length;
            int status;
            char* demangled = abi::__cxa_demangle(symbols[i], no_output_buffer, &length, &status);
//            std::cerr << i << " " << status << " ";
            std::cerr << (status == 0 ? demangled : info.dli_sname) << "\n";
            free(demangled);
        } else {
            std::cerr << i << " " << callstack[i] << "\n";
        }
    }
    free(symbols);
}

} // namespace Ariadne
