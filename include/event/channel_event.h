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


	// fragments
	int num;
	unsigned short daughter_charge[4];
	unsigned short daughter_mass[4];
	double daughter_energy[4];
	double daughter_time[4];
	double daughter_px[4];
	double daughter_py[4];
	double daughter_pz[4];
	double daughter_r[4];
	double daughter_theta[4];
	double daughter_phi[4];
	// beam
	double beam_energy;
	double beam_time;
	double beam_px;
	double beam_py;
	double beam_pz;
	double beam_r;
	double beam_theta;
	double beam_phi;
	// parent
	unsigned short parent_charge;
	unsigned short parent_mass;
	double parent_energy;
	double parent_px;
	double parent_py;
	double parent_pz;
	double parent_r;
	double parent_theta;
	double parent_phi;
	// recoil
	unsigned short recoil;
	unsigned short recoil_charge;
	unsigned short recoil_mass;
	double recoil_energy;
	double recoil_time;
	double recoil_px;
	double recoil_py;
	double recoil_pz;
	double recoil_r;
	double recoil_theta;
	double recoil_phi;
	// other information
	long long entry;
	int taf_index;
	int status[4];
};


}	// namespace ribll

#endif // __CHANNEL_EVENT_H__