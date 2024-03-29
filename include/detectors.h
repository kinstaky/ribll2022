#ifndef __DETECTORS_H__
#define __DETECTORS_H__

#include <memory>

#include "include/detector/adssd.h"
#include "include/detector/csi.h"
#include "include/detector/dssd.h"
#include "include/detector/ppac.h"
#include "include/detector/ssd.h"
#include "include/detector/tof.h"
#include "include/detector/vme_trigger.h"
#include "include/detector/t0d1.h"
#include "include/detector/t0d2.h"
#include "include/detector/t0d3.h"
#include "include/detector/tafd.h"
#include "include/detector/tabd.h"

namespace ribll {

/// @brief create detector by name
/// @param[in] name detector name
/// @param[in] run run number
/// @param[in] tag trigger tag
/// @returns pointer to detector if successful, nullptr otherwise
///
std::shared_ptr<Detector> CreateDetector(
	const std::string &name,
	unsigned int run,
	const std::string &tag
);


/// @brief create dssd detector by name
/// @param[in] name dssd name
/// @param[in] run run number
/// @param[in] tag trigger tag
/// @returns pointer to dssd if successful, nullptr otherwise
///
std::shared_ptr<Dssd> CreateDssd(
	const std::string &name,
	unsigned int run,
	const std::string &tag
);


/// @brief merge adssd trigger
/// @param[in] trigger_tag tag of XIA trigger
/// @param[in] run run number
/// @returns 0 if success, -1 otherwise
///
int MergeAdssdTrigger(const std::string &trigger_tag, unsigned int run);


/// @brief match events without trigger
/// @param[in] run run number
/// @returns 0 if success, -1 otherwise
///
int MatchWithoutTrigger(const std::string &detector, unsigned int run);


}

#endif 		// __DETECTORS_H__