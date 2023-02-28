#ifndef __PPAC_EVENT_H__
#define __PPAC_EVENT_H__

#include "include/defs.h"
#include "include/event/event.h"

namespace ribll {

class PpacMapEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	unsigned short index;
	unsigned short side;
	double time;
};


class PpacFundamentalEvent : public Event{
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
		return hit != 0;
	}


	/// @brief make the event invalid
	///
	virtual inline void Nullify() {
		hit = 0;
	}


	unsigned int flag;
	unsigned short hit;
	unsigned short x_hit;
	unsigned short y_hit;
	double x[ppac_num];
	double y[ppac_num];
};




}

#endif 		// __PPAC_EVENT_H__