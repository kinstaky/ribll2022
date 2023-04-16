#ifndef __PARTICLE_TYPE_EVENT_H__
#define __PARTICLE_TYPE_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class ParticleTypeEvent : public Event {
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


	unsigned short num;
	unsigned short charge[4];
	unsigned short mass[4];
	short layer[4];
};

}		// namespace ribll

#endif // __PARTICLE_TYPE_EVENT_H__