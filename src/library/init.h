/*
 * init.h
 *
 *  Created on: 18.12.2008
 *      Author: schwarz
 */

#ifndef INIT_H_
#define INIT_H_

#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <streambuf>
#include <cstring>
#include "TimeProf.h"

// Had to be additionally defined,
// because "char" otherwise will not be printed as a number.
//overloaded input operator
template<class T>
std::istream& operator>>(std::istream& is, std::vector<std::vector<T> >& tVec) {
  std::istringstream linestream;
  std::istringstream tokenstream;
  std::string line;
  std::string token;
  size_t pos, pcount, i;
  std::vector<T> result;

  if ((char) is.get() != '(')
    exit(1);
  std::getline(is, line);
  if (*line.rbegin() != ')')
    exit(1);
  linestream.str(line);
  tVec.clear();

  // if empty return to sender
	if ((line.length() > 0) && (line[0] == ')'))
		return is;

  while (!linestream.eof()) {
    std::getline(linestream, token);

#ifdef __DEBUG__
		std::cout << "token: " << token << std::endl;
#endif

    pos = pcount = 0;
    for (i = 0; i < token.size(); ++i) {
      if (token[i] == '(')
        ++pcount;
      if ((pcount == 0) && ((token[i] == ')') || (token[i] == ','))) {
        pos = i;
        break;
      }
      if (token[i] == ')')
        --pcount;
    }
    if (pos == 0) {
      std::cerr << "Error in file init.h, in function:" << std::endl;
			std::cerr << "std::istream& operator>>(std::istream& is, std::vector<std::vector<T> >& tVec)" << std::endl;
			std::cout << token << std::endl;
			throw;
		}

    tokenstream.clear();
    tokenstream.str(token.substr(0, pos));
    tokenstream >> result;
    tVec.push_back(result);

    linestream.clear();
    if (pos + 1 != token.size()) {
      linestream.str(token.substr(pos + 1, token.size() - pos - 1));
    } else {
      break;
    }

  }

  return is;
}

// Had to be additionally defined,
// because "char" otherwise will not be printed as a number.
//overloaded input operator
template<class T>
std::istream& operator>>(std::istream& is, std::vector<T>& tVec) {
  std::istringstream linestream;
  std::istringstream tokenstream;
  std::string line;
  std::string token;
  size_t pos, pcount, i;
  T result;
  double pre;

	// get first char
  if ((char) is.get() != '(')
    exit(1);

	// get line (without first char)
  std::getline(is, line);
  if (*line.rbegin() != ')')
    exit(1);
  linestream.str(line);
  tVec.clear();

  // if empty return to sender
	if ((line.length() > 0) && (line[0] == ')'))
		return is;

  while (!linestream.eof()) {
    std::getline(linestream, token);

#ifdef __DEBUG__
		std::cout << "token: " << token << std::endl;
#endif

    pos = pcount = 0;
    for (i = 0; i < token.size(); ++i) {
      if (token[i] == '(')
        ++pcount;
      if ((pcount == 0) && ((token[i] == ')') || (token[i] == ','))) {
        pos = i;
        break;
      }
      if (token[i] == ')')
        --pcount;
    }
    if (pos == 0) {
      std::cerr << "Error in file init.h, in function:" << std::endl;
			std::cerr << "std::istream& operator>>(std::istream& is, std::vector<T>& tVec)" << std::endl;
			throw;
		}

    tokenstream.clear();
    tokenstream.str(token.substr(0, pos));
    if (tokenstream.str() == "MAX") {
      tVec.push_back(std::numeric_limits<T>::max());
    } else {
      tokenstream >> pre;
      result = (T) pre;
      tVec.push_back(result);
    }

    linestream.clear();
    if (pos + 1 != token.size()) {
      linestream.str(token.substr(pos + 1, token.size() - pos - 1));
    } else {
      break;
    }

  }

  return is;
}

// Had to be additionally defined,
// because "char" otherwise will not be printed as a number.
template<class T>
std::ostream & operator<<(std::ostream &os, std::vector<std::vector<T> > &vec) {
  typename std::vector<std::vector<T> >::iterator it;

  os << "(";

  it = vec.begin();
  while (it != vec.end()) {
    os << *it;
    ++it;
    if (it != vec.end())
      os << ",";
  }

  os << ")";

  return os;
}

// Had to be additionally defined,
// because "char" otherwise will not be printed as a number.
template<class T>
std::ostream & operator<<(std::ostream &os, std::vector<T> &vec) {
  typename std::vector<T>::iterator it;

  os << "(";

  it = vec.begin();
  while (it != vec.end()) {
    if (*it == std::numeric_limits<T>::max())
      os << "MAX";
    else
      os << (double) *it;
    ++it;
    if (it != vec.end())
      os << ",";
  }

  os << ")";

  return os;
}

namespace rjungle {
  // this def.s make magic access fast
#ifdef __SPARSE_DATA__
  static const int magicAtMsk[4] = {192,48,12,3};
  static const int magicAtOfs[4] = {6,4,2,0};
  inline char magicAt(char* vec, size_t j) {
    return (vec[j / 4] & rjungle::magicAtMsk[j % 4]) >> rjungle::magicAtOfs[j % 4];
  }
#else
  template<class T>
  inline T magicAt(T* vec, uli_t j) {
    return vec[j];
  }
#endif

  static const char strLom[4][8] = {
    "undef", "nominal", "ordinal", "numeric" };
}

#endif /* INIT_H_ */
