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
	double be_channel[3];
	double he_channel[3];
	double ssd_channel[3];
	// position
	double be_x[3], be_y[3];
	double he_x[3], he_y[3];
	// recoil position
	double d_x, d_y;
	// target position
	double tx, ty;
	// PPAC flag
	unsigned short ppac_xflag, ppac_yflag;
	// PPAC position
	double ppac_x[3], ppac_y[3];
	// dssd normalize result information
	// dssd hit
	int be_x_hit[3], be_y_hit[3];
	int he_x_hit[3], he_y_hit[3];
	// dssd channel
	double be_x_channel[3][2], be_y_channel[3][2];
	double he_x_channel[3][2], he_y_channel[3][2];
	double d_x_channel, d_y_channel;
	// time
	double be_x_time[3][2], be_y_time[3][2];
	double he_x_time[3][2], he_y_time[3][2];
	double d_x_time, d_y_time;
	// strips
	unsigned int be_x_strip[3][2], be_y_strip[3][2];
	unsigned int he_x_strip[3][2], he_y_strip[3][2];
	unsigned int d_x_strip, d_y_strip;
	// state
	double q;
	// other information
	bool hole[2];
	// run
	int run;
};

} // ribll

#endif // __THREEBODY_INFO_EVENT_H__