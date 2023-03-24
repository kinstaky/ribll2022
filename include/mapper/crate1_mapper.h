#ifndef __CRATE1_MAPPER_H__
#define __CRATE1_MAPPER_H__

#include "include/mapper/xia_mapper.h"

namespace ribll {

class Crate1Mapper : public XiaMapper {
public:

	/// @brief constructor
	///
	/// @param[in] run run number
	///
	Crate1Mapper(unsigned int run);


	/// @brief default destructor
	///
	virtual ~Crate1Mapper() = default;


	/// @brief mapping function
	///
	/// @returns 0 for success, -1 otherwise
	///
	virtual int Map();
};

}		// namespace ribll

#endif		// __CRATE1_MAPPER_H__