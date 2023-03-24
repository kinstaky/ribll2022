#ifndef __CRATE2_MAPPER_H__
#define __CRATE2_MAPPER_H__

#include "include/mapper/xia_mapper.h"

namespace ribll {

class Crate2Mapper : public XiaMapper {
public:

	/// @brief constructor
	///
	/// @param[in] run run number
	///
	Crate2Mapper(unsigned int run);


	/// @brief default destructor
	///
	virtual ~Crate2Mapper() = default;


	/// @brief mapping function
	///
	/// @returns 0 for success, -1 otherwise
	///
	virtual int Map();
};

}		// namespace ribll

#endif		// __CRATE2_MAPPER_H__