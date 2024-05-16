#ifndef __PARTICLE_EVENT_H__
#define __PARTICLE_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class ParticleEvent : public Event {
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
	unsigned short charge[8];
	unsigned short mass[8];
	double energy[8];
	double time[8];
	double x[8];
	double y[8];
	double z[8];
	double px[8];
	double py[8];
	double pz[8];
	int status[8];
	int index[8];
};

}	// namespace ribll

#endif	// __PARTICLE_EVENT_H__