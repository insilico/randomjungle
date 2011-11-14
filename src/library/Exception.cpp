#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "Exception.h"

Exception::Exception(std::string str): str(str), status(0) {
  std::cerr << "ERROR:" << std::endl << str << std::endl;
}

const char* Exception::what() const throw() {
  std::string strret("");

  if (status != 0) {
    char result[100];
    sprintf(result, "sta:%d ", status);
    strret.append(result);
  }

  return strret.c_str();
}
