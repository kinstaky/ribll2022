#ifndef __CSI_EVENT_H__
#define __CSI_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class CsiMapEvent : public Event {
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


	unsigned short index;
	bool cfd_flag;
	double time;
	double energy;
	long long decode_entry;
};


class CircularCsiFundamentalEvent : public Event {
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


	bool match;
	unsigned short cfd_flag;
	double time[12];
	double energy[12];
	long long decode_entry[12];
};


class SquareCsiFundamentalEvent : public Event {
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


	bool match;
	unsigned short cfd_flag;
	double time[4];
	double energy[4];
	long long decode_entry[4];
};

}

#endif 		// __CSI_EVENT_H__