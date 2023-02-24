#ifndef __PPAC_EVENT_H__
#define __PPAC_EVENT_H__

#include "include/defs.h"
#include "include/event/event.h"

namespace ribll {

struct PpacMapEvent : public Event {
	unsigned short index;
	unsigned short side;
	double time;


	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;
};


struct PpacFundamentalEvent : public Event{
	unsigned int flag;
	unsigned short hit;
	unsigned short x_hit;
	unsigned short y_hit;
	double x[ppac_num];
	double y[ppac_num];


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

#endif 		// __PPAC_EVENT_H__