#ifndef __PD_INFO_EVENT_H__
#define __PD_INFO_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class PdInfoEvent : public Event {
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
	// energy
	double tafd_energy;
	double taf_energy;
	double t0_energy;
	// channel
	double csi_channel;
	double c_channel[2];
	// position
	double c_x[2], c_y[2];
	// recoil position
	double d_x, d_y;
	// target position
	double tx, ty;
	// PPAC flag, 0-XPPAC, 1-VPPAC
	int ppac_flag;
	// PPAC x,y flag
	unsigned short ppac_xflag, ppac_yflag;
	// PPAC position
	double ppac_x[3], ppac_y[3];
	// PPAC number used for tracking
	int ppac_x_track, ppac_y_track;
	// dssd normalize result information
	// // dssd hit
	// int c_x_hit[2], c_y_hit[2];
	// // dssd channel
	// double c_x_channel[2][2], c_y_channel[2][2];
	// double d_x_channel, d_y_channel;
	// // time
	// double c_x_time[2][2], c_y_time[2][2];
	// double d_x_time, d_y_time;
	// // strips
	// unsigned int c_x_strip[2][2], c_y_strip[2][2];
	// unsigned int d_x_strip, d_y_strip;
	// state
	double q;
	double c15_kinetic;
	double c15_ex;
	// run
	int run;
	long long entry;
};

} // ribll

#endif // __PD_INFO_EVENT_H__