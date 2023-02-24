#ifndef __TOF_EVENT_H__
#define __TOF_EVENT_H__

#include "include/event/event.h"

namespace ribll {

struct TofMapEvent : public Event {
	double time;
	unsigned short index;


	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;
};

struct TofFundamentalEvent : public Event {
	double time[2];


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

#endif 		// __TOF_EVENT_H__