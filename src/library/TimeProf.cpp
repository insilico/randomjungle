#include "TimeProf.h"


TimeProf::TimeProf() {
}

TimeProf::~TimeProf() {
}


void TimeProf::start(std::string name) {
	std::map<std::string, Profiler>::iterator it;

	// search profiler in map
	it = profilers.find(name);

	// add a new if not there
	if (it == profilers.end())
		profilers[name] = Profiler();

	// start timer
	profilers[name].start();
}


void TimeProf::stop(std::string name) {
	std::map<std::string, Profiler>::iterator it;

	// search profiler in map
	it = profilers.find(name);

	// profiler has to be there
	assert(it != profilers.end());

	// stop timer
	profilers[name].stop();
}

void TimeProf::writeToFile(std::ostream &io) {
	std::map<std::string, Profiler>::iterator it;

	// get first profiler
	it = profilers.begin();

	// show all profilers
	while(it != profilers.end()) {

		// print name and total time that was consumed
		io << it->first << " " << it->second.totalTime << std::endl;

		++it;
	}
}

#ifdef HAVE_TIMEPROF

// init one global
TimeProf timeProf;

#endif

