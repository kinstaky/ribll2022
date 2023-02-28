#ifndef __DSSD_EVENT_H__
#define __DSSD_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class DssdMapEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	unsigned short side;
	unsigned short strip;
	double time;
	double energy;
};


class DssdFundamentalEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	/// @brief show whether the event is valid or not
	/// @returns true if valid, false if invalid
	///
	virtual inline bool Valid() const {
		return front_hit != 0 && back_hit != 0;
	}


	/// @brief make the event invalid
	///
	virtual inline void Nullify() {
		front_hit = back_hit = 0;
	}
	

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

#endif		// __DSSD_EVENT_H__