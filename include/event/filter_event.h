#ifndef __FILTER_EVENT_H__
#define __FILTER_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class FilterEvent : public Event{
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
	unsigned short pid_index[4];
	unsigned short merge_flag[4];
	unsigned short norm_front_index[4];
	unsigned short norm_back_index[4];
	unsigned short front_index[4];
	unsigned short back_index[4];
	unsigned short front_strip[4];
	unsigned short back_strip[4];
	double front_energy[4];
	double back_energy[4];
};

} // namespace ribll

#endif // __FILTER_EVENT_H__