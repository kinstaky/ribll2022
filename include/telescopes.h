#ifndef __TELECOPES_H__
#define __TELECOPES_H__

#include <memory>

#include "include/telescope/taf.h"

namespace ribll {

/// @brief create telescope by name
/// @param[in] name telescope name
/// @param[in] run run number
/// @param[in] tag trigger tag
/// @returns pointer to telescope if successful, nullptr otherwise
///
std::shared_ptr<Telescope> CreateTelescope(
	const std::string &name,
	unsigned int run,
	const std::string &tag
);

}		// namespace ribll

#endif		// __TELECOPES_H__