#ifndef __SSD_EVENT_H__
#define __SSD_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class SsdEvent : public Event {
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


	bool cfd_flag;
	double time;
	double energy;
};

}		// namespace ribll

#endif 		// __SSD_EVENT_H__