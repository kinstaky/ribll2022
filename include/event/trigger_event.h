#ifndef __TRIGGER_EVENT_H__
#define __TRIGGER_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class TriggerEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	/// @param[in] prefix prefix of variables for friend tree
	///
	virtual void SetupInput(
		TTree *tree,
		const std::string &prefix = ""
	) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	double time;
	long long timestamp;
	bool cfd_flag;
};

}

#endif		// __TRIGGER_EVENT_H__