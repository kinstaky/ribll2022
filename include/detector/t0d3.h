#ifndef __T0D3_H__
#define __T0D3_H__

#include "include/detector/dssd.h"

namespace ribll {

class T0d3 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	///
	T0d3(unsigned int run);


	/// @brief default destructor
	///
	virtual ~T0d3() = default;
};

}		// namespace ribll

#endif		// __T0D3_H__