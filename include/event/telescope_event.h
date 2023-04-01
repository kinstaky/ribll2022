#ifndef __TELECOPE_EVENT_H__
#define __TELECOPE_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class TelescopeEvent : public Event {
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


	unsigned short particle;
	unsigned short layer[4];
	unsigned short flag[4];
	double energy[4][8];
	double radius[4][8];
	double theta[4][8];
	double phi[4][8];
};


}		// namespace ribll

#endif		// __TELECOPE_EVENT_H__