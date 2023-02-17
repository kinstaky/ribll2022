#ifndef __EVENT_H__
#define __EVENT_H__

#include "include/defs.h"

#include <TTree.h>

namespace ribll {

struct Event {
	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) = 0;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) = 0;
};

// struct Event {
// 	long long timestamp;
// };

// struct SingleEvent : Event {
// 	unsigned short index;
// 	unsigned short side;
// 	unsigned short strip;
// 	double time;
// 	double energy;
// };

struct PPACEvent {
	long long timestamp;
	// flags and hits
	int flag;
	unsigned short hit;
	unsigned short x_hit;
	unsigned short y_hit;
	// time in nanoseconds
	double x[ppac_num];				// x1-x2
	double y[ppac_num];				// y1-y2
};

struct DSSDEvent {
	long long timestamp;
	unsigned short index;
	unsigned short front_hit;
	unsigned short back_hit;
	unsigned short front_strip[8];
	unsigned short back_strip[8];
	double front_time[8];
	double back_time[8];
	double front_energy[8];
	double back_energy[8];
};

}

#endif		// __EVENT_H__