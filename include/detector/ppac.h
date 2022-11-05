#ifndef __PPAC_H__
#define __PPAC_H__

#include <map>
#include <vector>

#include <TROOT.h>

namespace ribll {

class PPAC {
public:

	/// @brief constructor
	///
	/// @param[in] run
	/// 
	PPAC(int run);


	/// @brief correlate and build hit events
	///
	/// @returns 0 for success, -1 otherwise
	///
	int Correlate();

private:
	struct Event {
		unsigned short index;
		unsigned short side;
		double time;
		bool used;
	};

	int run_;
	Event event_;

	std::vector<Long64_t> trigger_time_;
	std::multimap<Long64_t, Event> events_;

	// int ReadTrigger();

	int ReadEvents();
};

}

#endif