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


	int num;
	short layer[8];
	unsigned short flag[8];
	unsigned short ssd_flag;
	unsigned short charge[8];
	unsigned short mass[8];
	double energy[8][3];
	double time[8][3];
	double ssd_energy[3];
	double ssd_time[3];
	double x[8][3];
	double y[8][3];
	double z[8][3];
	int status[8];
	int points[8];
	unsigned short dssd_flag[8][3];
	bool hole[8];
};

}	// namespace ribll

#endif		// __T0_EVENT_H__