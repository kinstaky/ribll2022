#ifndef __GENERATE_EVENT_H__
#define __GENERATE_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class GenerateEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pinter to input tree
	/// @param[in] prefix prefix of variables for friend tree
	///
	virtual void SetupInput(
		TTree *tree,
		const std::string &prefix = ""
	) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pinter to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	double beam_kinetic;
	double beam_excited_energy;
	double beam_theta;
	double beam_phi;
	double target_x;
	double target_y;
	double fragment_excited_energy;
	int fragment_state;
	double parent_kinetic;
	double parent_theta;
	double parent_phi;
	double recoil_kinetic;
	double recoil_theta;
	double recoil_phi;
	double rx;
	double ry;
	double rz;
	double rr;
	double fragment_kinetic[2];
	double fragment_theta[2];
	double fragment_phi[2];
	double fragment_x[2];
	double fragment_y[2];
	double fragment_z[2];
	double fragment_r[2];
};

}	// ribll2022

#endif	//__GENERATE_EVENT_H__