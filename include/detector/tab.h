#ifndef __TAB_H__
#define __TAB_H__

#include "include/detector/adssd.h"

namespace ribll {

class Tab : public Adssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] index index of taf, 0 to 5
	/// @param[in] tag trigger tag
	///
	Tab(unsigned int run, unsigned int index, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~Tab() = default;
};

}

#endif		// __TAB_H__