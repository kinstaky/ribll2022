#ifndef __TOF_EVENT_H__
#define __TOF_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class TofMapEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	double time;
	unsigned short index;
};


class TofFundamentalEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	double time[2];
};


}

#endif 		// __TOF_EVENT_H__