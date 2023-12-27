#ifndef __DETECT_EVENT_H__
#define __DETECT_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class DetectEvent : public Event {
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


	int taf_layer;
	double taf_lost_energy[2];
	double taf_energy[2];
	// particle stop layer, 0-T0D1, 1-T0D2, 2-T0D3
	int t0_layer[2];
	double t0_lost_energy[2][7];
	double t0_energy[2][7];
	double tafx;
	double tafy;
	double tafz;
	double tafr;
	double t0x[2][3];
	double t0y[2][3];
	double t0z[2][3];
	double t0r[2][3];
	double ppacx[3];
	double ppacy[3];
	// target x from PPAC
	double tx;
	// target y from PPAC
	double ty;
	int valid;
	double q;
};


}	// ribll

#endif // __DETECT_EVENT_H__