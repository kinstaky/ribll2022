#ifndef __DSSD_EVENT_H__
#define __DSSD_EVENT_H__

#include "include/event/event.h"

namespace ribll {

struct DssdMapEvent : public Event {
	unsigned short index;
	unsigned short side;
	unsigned short strip;
	double time;
	double energy;


	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;
};


struct DssdFundamentalEvent : public Event {
	unsigned short index;
	unsigned short front_hit;
	unsigned short back_hit;
	unsigned short front_strip[8];
	unsigned short back_strip[8];
	double front_time[8];
	double back_time[8];
	double front_energy[8];
	double back_energy[8];


	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;
};

}

#endif		// __DSSD_EVENT_H__