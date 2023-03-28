#ifndef __T0D2_H__
#define __T0D2_H__

#include "include/detector/dssd.h"

namespace ribll {

class T0d2 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	///
	T0d2(unsigned int run, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~T0d2() = default;
};

}		// namespace ribll

#endif		// __T0D2_H__