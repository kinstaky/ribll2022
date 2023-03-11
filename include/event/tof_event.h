#ifndef __TOF_EVENT_H__
#define __TOF_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class TofMapEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	double time;
	unsigned short index;
	bool cfd_flag;
};


class TofFundamentalEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	///
	virtual void SetupInput(TTree *tree) override;


	/// @brief setup branches of output tree
	/// @param[out] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	/// @brief show whether the event is valid or not
	/// @returns true if valid, false if invalid
	///
	virtual inline bool Valid() const {
		return time[0] > -9e4 && time[1] > -9e4;
	}


	/// @brief make the event invalid
	///
	virtual inline void Nullify() {
		time[0] = time[1] = -1e5;
	}


	double time[2];
	unsigned short cfd_flag;
};


}

#endif 		// __TOF_EVENT_H__