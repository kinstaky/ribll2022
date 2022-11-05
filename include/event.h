#ifndef __EVENT_H__
#define __EVENT_H__

#include "include/defs.h"

namespace ribll {

struct Event {
	long long timestamp;
};

struct SingleEvent : Event {
	unsigned short index;
	unsigned short side;
	unsigned short strip;
	double time;
	double energy;
};

struct PPACEvent : Event {
	// flags and hits
	int flag;
	unsigned short hit;
	unsigned short x_hit;
	unsigned short y_hit;
	// time in nanoseconds
	double x[ppac_num];				// x1-x2
	double y[ppac_num];				// y1-y2
};

struct DSSDEvent : Event {
	unsigned short index;
	unsigned short x_hit;
	unsigned short y_hit;
	unsigned short x_strip[8];
	unsigned short y_strip[8];
	double x_time[8];
	double y_time[8];
	double x_energy[8];
	double y_energy[8];
};

}

#endif		// __EVENT_H__