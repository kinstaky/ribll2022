#ifndef __T0_EVENT_H__
#define __T0_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class T0Event : public Event {
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
	unsigned short layer[4];
	unsigned short flag[4];
	unsigned short ssd_flag;
	double energy[4][4];
	double ssd_energy[3];
	double x[4][4];
	double y[4][4];
	double z[4][4];
};

}	// namespace ribll

#endif		// __T0_EVENT_H__