/**
 * @author Roman Pahl (original) and Daniel F Schwarz (modified heavily)
 */

#ifndef IOLINES_H_
#define IOLINES_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>
#include "Exception.h"

/**
 * @brief A class for iterating through tokens from a c-string
 */
class StringIterator {
public:
  StringIterator(char delimiters[] = " \n\r") :
    strTmp(NULL), myLastToken(NULL) {
    strcpy(this->delimiters, delimiters);
  }

  StringIterator() :
    strTmp(NULL), myLastToken(NULL) {
  }

  virtual ~StringIterator() {
    // free string memory
    if (strTmp)
      delete strTmp;
  }

  inline void setDelimiters(char delimiters[] = " \n\r") {
    strcpy(this->delimiters, delimiters);
  }

  /**
   * @brief Initializes the token
   * @param theString the string to extract from
   */
  char* init(char *theString) {
    token = mystrtok(theString, delimiters);

    return token;
  }

  /**
   * @brief Initializes the token
   * @param theString the string to extract from
   */
  char* init(std::string theCppString) {
    // free string memory
    if (strTmp)
      delete strTmp;

    // copy string
    strTmp = new char[strlen(theCppString.c_str()) + 1];
    strcpy(strTmp, theCppString.c_str());

    // get token
    token = mystrtok(strTmp, delimiters);

    return token;
  }

  //! @return true if token is empty, else false
  inline bool empty() const {
    return token == NULL;
  }

  //! Jump to next token
  char* next() {
    token = mystrtok(NULL, delimiters);
    return token;
  }

  //! Jump to next token
  char* next(char* delims) {
    token = mystrtok(NULL, delims);
    return token;
  }

  //! Get next string
  char *mystrtok(char *s1, const char *delimit) {
    char *tmp;

    /* Skip leading delimiters if new string. */
    if (s1 == NULL) {
      s1 = myLastToken;
      if (s1 == NULL) /* End of story? */
        return NULL;
    } else {
      s1 += strspn(s1, delimit);
    }

    /* Find end of segment */
    tmp = strpbrk(s1, delimit);
    if (tmp) {
      /* Found another delimiter, split string and save state. */
      this->lastDelim = *tmp;
      *tmp = '\0';
      myLastToken = tmp + 1;
    } else {
      /* Last segment, remember that. */
      myLastToken = NULL;
    }

    return s1;
  }

  //! The token
  char *token;

  //! temp string
  char *strTmp;

  //! last delimiter
  char lastDelim;

  //! last token
  char *myLastToken;

private:
  // Only needed if delimiters are chosen at runtime
  char delimiters[10]; // modify here if you want more than 10 delimiters
};

#endif /* IOLINES_H_ */
