// Copyright 2008, 2009 Martin C. Frith

// Convert things (e.g. numbers) to and from strings.  These are
// replacements for boost::lexical_cast: by avoiding dependency on
// boost, we can distribute code more easily.

#ifndef STRINGIFY_HH
#define STRINGIFY_HH
#include <sstream>
#include <string>
#include <stdexcept>
#include <cassert>

namespace cbrc{

template<typename T>
std::string stringify( const T& x ){
  std::ostringstream oss;
  oss << x;
  assert(oss);
  return oss.str();
}

template<typename T>
void unstringify( T& x, const std::string& s ){
  std::istringstream iss(s);
  if( !(iss >> x) || !(iss >> std::ws).eof() ){
    throw std::runtime_error( "can't interpret: " + s );
  }
}

template<typename T>
void unstringifySize( T& x, const std::string& s ){
  std::istringstream iss(s);
  if( !(iss >> x) ){
    throw std::runtime_error( "can't interpret: " + s );
  }

  std::string suffix;
  if( iss >> suffix ){
    int i;
    /**/ if( suffix == "K" ) i = 1;  // "KibiBytes"
    else if( suffix == "M" ) i = 2;  // "MebiBytes"
    else if( suffix == "G" ) i = 3;  // "GibiBytes"
    else if( suffix == "T" ) i = 4;  // "TebiBytes"
    else if( suffix == "P" ) i = 5;  // "PebiBytes"
    else throw std::runtime_error( "can't interpret: " + s );

    while( i-- ){
      if( (x * 1024) / 1024 != x ){
	throw std::runtime_error( "can't interpret (too big): " + s );
      }
      x *= 1024;
    }
  }

  if( !(iss >> std::ws).eof() ){
    throw std::runtime_error( "can't interpret: " + s );
  }
}

}

#endif
