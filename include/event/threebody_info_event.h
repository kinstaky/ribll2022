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
	// TAF flag
	//   0 - normal event with dE-E PID, particle stopped in CsI
	//   1 - particle stopped in ADSSD
	//   2 - particle stopped in CsI, but cut by CsI threshold
	int taf_flag;
	// energy
	double t0_energy[2];
	double tafd_energy;
	double csi_energy;
	double taf_energy;
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
	// target position from XIA PPAC
	double xptx, xpty;
	// target position from VME PPAC
	double vptx, vpty;
	// PPAC flag, XPPAC-bit_0, VPPAC-bit_1
	int ppac_flag;
	// XIA PPAC x,y flag
	unsigned short xppac_xflag, xppac_yflag;
	// VME PPAC x,y flag
	unsigned short vppac_xflag, vppac_yflag;
	// XIA PPAC position
	double xppac_x[3], xppac_y[3];
	// VME PPAC position
	double vppac_x[3], vppac_y[3];
	// XIA PPAC track number, index 0 for x, 1 for y
	int xppac_track[2];
	// VME PPAC track number, index 0 for x, 1 for y
	int vppac_track[2];
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
	double c14_kinetic;
	// Q from XIA PPAC reaction point or VME PPAC reaction point
	double q, vq;
	// other information
	// bind flag, show bind strips side in bits
	// bit 1-D1F, bit 2-D1B, bit 3-D2F, bit 4-D2B
	// e.g. 0x9 means T0D1 is F1B2 and T0D2 is F2B1, 0x0 means no bind strips
	int bind;
	// T0D2 hole flag
	bool hole[2];
	// run
	int run;
	long long entry;
};

} // ribll

#endif // __THREEBODY_INFO_EVENT_H__