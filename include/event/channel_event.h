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
	unsigned short charge[4];
	unsigned short mass[4];
	double energy[4];
	double px[4];
	double py[4];
	double pz[4];
	double r[4];
	double theta[4];
	double phi[4];
};


}	// namespace ribll

#endif // __CHANNEL_EVENT_H__