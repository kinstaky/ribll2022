#ifndef __TA_EVENT_H__
#define __TA_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class TaEvent : public Event {
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


	int num;
	unsigned short flag[4];
	double energy[4][4];
	double time[4][4];
	double radius[4][4];
	double theta[4][4];
	double phi[4][4];
};


}		// namespace ribll

#endif		// __TA_EVENT_H__