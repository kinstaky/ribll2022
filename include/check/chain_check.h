#ifndef __CHAIN_CHECKER__
#define __CHAIN_CHECKER__

#include <string>

namespace ribll {


/// @brief check whether the key is valid, i.e. can be found in dependent map
/// @param[in] key key to check
/// @returns true if the key is valid, false otherwise
///
bool NodeCheck(const std::string &key);


/// @brief check whether the process is outdated
/// @param[in] run run number
/// @param[in] key process key, used in dependency map
/// @returns true if the key is up to date, false if outdated
///
bool ChainCheck(unsigned int run, const std::string &key);


}	// namespace ribll

#endif		// __CHAIN_CHECKER__