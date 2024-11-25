#ifndef __CHANNEL_V2_EVENT_H__
#define __CHANNEL_V2_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class ChannelV2Event : public Event {
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

	// valid flags
	// valid bits: 0: T0, 1: TAF, 2: PPAC
	int valid;
	bool ppac_valid;
	int hole;
	bool t0_valid;
	bool tafd_edge;
	bool tafcsi_valid;	
	// PPAC information
	unsigned short ppac_xflag, ppac_yflag;
	unsigned short ppac_xnum, ppac_ynum;
	double ppac_x[3], ppac_y[3];
	// target position
	double tx, ty;
	// beam particle
	unsigned short beam_charge, beam_mass;
	double beam_energy, beam_kinetic, beam_momentum;
	// parent particle
	unsigned short parent_charge, parent_mass;
	double parent_momentum;
	// recoil particle
	int taf_index, csi_index;
	unsigned short recoil_charge, recoil_mass;
	double recoil_x, recoil_y, recoil_z;
	double recoil_energy, recoil_kinetic, recoil_momentum;
	// fragments
	int fragment_num;
	unsigned short fragment_charge[2], fragment_mass[2];
	double fragment_x[2], fragment_y[2], fragment_z[2];
	double fragment_energy[2], fragment_kinetic[2], fragment_momentum[2];
	// other information
	short run;
	long long entry;
	// extra TAF information
	int tafd_front_strip;
	double tafd_energy;
};

}	// ribll

#endif // __CHANNEL_V2_EVENT_H__