#ifndef __THREEBODY_INFO_EVENT_H__
#define __THREEBODY_INFO_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class ThreeBodyInfoEvent : public Event {
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


	// index
	int csi_index;
	// layers
	int layer[2];
	// energy
	double tafd_energy;
	double t0_energy[2];
	// channel
	double csi_channel;
	double d1_channel[2];
	double d2_channel[2];
	double d3_channel[2];
	double ssd_channel[3];
	// position
	double d1x[2];
	double d1y[2];
	double d2x[2];
	double d2y[2];
	double d3x[2];
	double d3y[2];
	// recoil position
	double rx, ry;
	// target position
	double tx, ty;
	// PPAC flag
	unsigned short ppac_xflag, ppac_yflag;
	// PPAC position
	double ppac_x[3], ppac_y[3];
	// run
	int run;
};

} // ribll

#endif // __THREEBODY_INFO_EVENT_H__