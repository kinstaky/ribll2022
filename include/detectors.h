#ifndef __DETECTORS_H__
#define __DETECTORS_H__

#include <memory>

#include "include/detector/tof.h"
#include "include/detector/vme_trigger.h"
#include "include/detector/dssd.h"
#include "include/detector/adssd.h"
#include "include/detector/csi.h"

namespace ribll {

/// @brief create detector by name
/// @param[in] name detector name 
/// @param[in] run run number 
/// @returns pointer to detector if successful, nullptr otherwise
///
std::shared_ptr<Detector> CreateDetector(
	const std::string &name,
	unsigned int run
);

}

#endif 		// __DETECTORS_H__