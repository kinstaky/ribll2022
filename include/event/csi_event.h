#ifndef __CSI_EVENT_H__
#define __CSI_EVENT_H__

#include "include/event/event.h"

namespace ribll {

struct CsiMapEvent : public Event {
	unsigned short index;
	double time;
	double energy;


	/// @brief setup branches of input tree
	/// @param[in] tree pinter to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pinter to output tree
	///
	virtual void SetupOutput(TTree *tree) override;
};


struct CircularCsiFundamentalEvent : public Event {
	bool match;
	double time[12];
	double energy[12];


	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param tree pointer to output tree
	virtual void SetupOutput(TTree *tree) override;
};


struct SquareCsiFundamentalEvent : public Event {
	bool match;
	double time[4];
	double energy[4];


	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param tree pointer to output tree
	virtual void SetupOutput(TTree *tree) override;
};

}

#endif 		// __CSI_EVENT_H__