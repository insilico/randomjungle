#ifndef EXCEPTION_H_
#define EXCEPTION_H_

/*
 * Includes
 */

//#include <stdexcept>

#include <string>

/*
 * Class declaration
 */

class Exception //: public std::exception
{
public:

  Exception(std::string str);
	Exception(int status, std::string str): str(str), status(status) {};

	virtual ~Exception() {};

	virtual const char* what() const throw();

	std::string str;
	int status;

};

#endif /*EXCEPTION_H_*/
