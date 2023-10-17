// Copyright 2008, 2009, 2010 Martin C. Frith

#ifndef IO_HH
#define IO_HH

#include <string>
#include <stdexcept>
#include <fstream>
#include <cassert>

namespace cbrc{

template<typename T>  // T should be a vector-iterator or a pointer
void memoryToStream( T beg, T end, std::ostream& s ){
  assert( beg < end );
  enum { CHUNK_SIZE = 1073741824 };  // need to do big writes in chunks: why?
  const char * b = (const char*)&(*beg);
  const char * e = (const char*)&(*end);
  while( e - b > CHUNK_SIZE ){
    s.write( b, CHUNK_SIZE );
    b += CHUNK_SIZE;
  }
  s.write( b, e - b );
}

template<typename T>  // T should be a vector-iterator or a pointer
void memoryToBinaryFile( T beg, T end, const std::string& fileName ){
  if( beg == end ) return;
  std::ofstream file( fileName.c_str(), std::ios::binary );
  if( !file ) throw std::runtime_error( "can't open file: " + fileName );
  memoryToStream( beg, end, file );
  file.close();
  if( !file ) throw std::runtime_error( "can't write file: " + fileName );
}

}

#endif
