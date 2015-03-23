# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

project(tdigest)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
include(cmake/dist.cmake)
include(lua)

if(MSVC)
    # Predefined Macros: http://msdn.microsoft.com/en-us/library/b0084kay.aspx
    # Compiler options: http://msdn.microsoft.com/en-us/library/fwkeyyhe.aspx
    # set a high warning level and treat them as errors
    set(CMAKE_CXX_FLAGS           "/W3 /WX")

    # debug multi threaded dll runtime, complete debugging info, runtime error checking
    set(CMAKE_CXX_FLAGS_DEBUG     "/MDd /Zi /RTC1")

    # multi threaded dll runtime, optimize for speed, auto inlining
    set(CMAKE_CXX_FLAGS_RELEASE   "/MD /O2 /Ob2 /DNDEBUG")
else()
    # Predefined Macros: clang|gcc -dM -E -x c /dev/null
    # Compiler options: http://gcc.gnu.org/onlinedocs/gcc/Invoking-GCC.html#Invoking-GCC
    set(CMAKE_CXX_FLAGS   "-std=c++0x -pedantic -Werror -Wall -Wextra")
    if (NOT WIN32)
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fPIC")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG     "-g")
    set(CMAKE_CXX_FLAGS_RELEASE   "-O2 -DNDEBUG")
endif()

include_directories(include)
add_definitions ( -DDIST_VERSION="${DIST_VERSION}" )
install_lua_module(tdigest src/TDigest.cpp src/lua.cpp src/TDigest.def )

add_lua_test(tests/test.lua)
