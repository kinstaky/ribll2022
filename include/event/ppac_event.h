#ifndef __PPAC_EVENT_H__
#define __PPAC_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class PpacMapEvent : public Event {
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


	unsigned short index;
	unsigned short side;
	bool cfd_flag;
	double time;
};


class PpacFundamentalEvent : public Event{
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


	/// @brief show whether the event is valid or not
	/// @returns true if valid, false if invalid
	///
	virtual inline bool Valid() const {
		return hit != 0;
	}


	/// @brief make the event invalid
	///
	virtual inline void Nullify() {
		hit = 0;
		flag = 0;
	}


	unsigned int flag;
	unsigned short cfd_flag;
	unsigned short hit;
	unsigned short x_hit;
	unsigned short y_hit;
	double x1[3];
	double x2[3];
	double y1[3];
	double y2[3];
	double anode[3];
};


class PpacMergeEvent : public Event {
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


	unsigned short xflag;
	unsigned short yflag;
	double x[3];
	double y[3];
	double time[3];
};




}

#endif 		// __PPAC_EVENT_H__