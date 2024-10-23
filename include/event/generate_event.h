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


	double depth;
	double beam_kinetic_in_target;
	double beam_kinetic_before_target;
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
	double recoil_kinetic_after_target;
	double recoil_kinetic_in_target;
	double recoil_theta;
	double recoil_phi;
	double rx;
	double ry;
	double rz;
	double rr;
	double fragment_kinetic_after_target[2];
	double fragment_kinetic_in_target[2];
	double fragment_theta[2];
	double fragment_phi[2];
	double fragment_x[2];
	double fragment_y[2];
	double fragment_z[2];
	double fragment_r[2];
	// other information
	double elastic_angle;
	double breakup_angle;
	double parent_recoil_angle;
	double fragment_phi_center;
	double fragment_fragment_angle;
	double angle_theta_star;
	double angle_psi;
	double angle_chi;
	// velocity information
	double recoil_vx, recoil_vy, recoil_vz;
	double parent_vx, parent_vy, parent_vz;
	double fragment_vx[2], fragment_vy[2], fragment_vz[2];
};

}	// ribll2022

#endif	//__GENERATE_EVENT_H__