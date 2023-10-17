// Copyright 2010, 2011 Martin C. Frith

// Functions for mapping files into memory.

#ifndef FILE_MAP_HH
#define FILE_MAP_HH

#include <stddef.h>  // size_t
#include <string>

namespace cbrc{

// Maps a file into memory, read-only, and returns a pointer to the
// start of the mapping.  If it fails, it throws a runtime_error.  If
// bytes is zero, it does nothing and returns 0.
void* openFileMap( const std::string& fileName, size_t bytes );

// Releases a file mapping.  If bytes is zero, it does nothing.  If it
// fails, it throws a runtime_error.
void closeFileMap( void* begin, size_t bytes );

}

#endif
