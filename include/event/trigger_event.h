#ifndef __TRIGGER_EVENT_H__
#define __TRIGGER_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class TriggerEvent : public Event {
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
};

}

#endif		// __TRIGGER_EVENT_H__