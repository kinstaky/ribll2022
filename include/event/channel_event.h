#ifndef __CHANNEL_EVENT_H__
#define __CHANNEL_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class ChannelEvent : public Event {
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


	unsigned short num;
	unsigned short charge[8];
	unsigned short mass[8];
	double energy[8];
	double px[8];
	double py[8];
	double pz[8];
	double r[8];
	double theta[8];
	double phi[8];
};


}	// namespace ribll

#endif // __CHANNEL_EVENT_H__