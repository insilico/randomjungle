/* 
 * Copyright (C) 2008-2010  Daniel F. Schwarz
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TIMEPROF_H_
#define TIMEPROF_H_

// define class if time profiling is wanted

#include <map>
#include <ctime>
#include <cassert>
#include <iostream>

#include "Profiler.h"
#include "config.h"

class TimeProf {
	public:

		/** 
		 * \brief Constructor
		 */
		TimeProf();

		/** 
		 * \brief Destructor
		 */
		virtual ~TimeProf();

		/** 
		 * \brief Start the stop watch
		 */
		void start(std::string name);

		/** 
		 * \brief Stop the watch and add time difference to total time
		 */
		void stop(std::string name);

		//! Streams for output profiled informations 
		void writeToFile(std::ostream &io);

		//! Map of all time trigger and profiler
		std::map<std::string, Profiler> profilers;
};

#ifdef HAVE_TIMEPROF

extern TimeProf timeProf;

#define TIMEPROF_START(x) timeProf.start(x);
#define TIMEPROF_STOP(x) timeProf.stop(x);

#else /*HAVE_TIMEPROF*/

// define empty macros
#define TIMEPROF_START(x) ;
#define TIMEPROF_STOP(x) ; 

#endif /*HAVE_TIMEPROF*/

#endif /*TIMEPROF_H_*/
