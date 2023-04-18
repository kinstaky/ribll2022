#ifndef __GLOBAL_DEF_H__
#define __GLOBAL_DEF_H__

#include <vector>

namespace ribll {

//-----------------------------------------------------------------------------
//									data path
//-----------------------------------------------------------------------------
// decode data
const char* const kCrate0Path = "/data2/ribll2022/decode/crate1/";
const char* const kCrate0FileName = "c1data";

const char* const kCrate1Path = "/data2/ribll2022/decode/crate2/";
const char* const kCrate1FileName = "c2data";

const char* const kCrate2Path = "/data2/ribll2022/decode/crate3/";
const char* const kCrate2FileName = "c3data";

const char* const kCrate3Path = "/data2/ribll2022/decode/vme/";
const char* const kCrate3FileName = "vmedata";

// generate data
const char* const kGenerateDataPath = "/data2/ribll2022/";
const char* const kMappingDir = "map/";
const char* const kDecodeDir = "decode/";
const char* const kAlignDir = "align/";
const char* const kFundamentalDir = "fundamental/";
const char* const kBeamDir = "beam/";
const char* const kNormalizeDir = "normalize/";
const char* const kMergeDir = "merge/";
const char* const kTelescopeDir = "telescope/";
const char* const kCalibrationDir = "calibration/";
const char* const kParticleIdentifyDir = "pid/";
const char* const kParticleDir = "particle/";

const char* const kCheckDir = "check/";
const char* const kShowDir = "show/";
const char* const kEnergyCalculateDir = "energy_calculate/";


//-----------------------------------------------------------------------------
//								xia vme alignment
//-----------------------------------------------------------------------------
const int kScalerIndex = 4;
const int kScalerPeriod = 200;

//-----------------------------------------------------------------------------
//								vme module config
//-----------------------------------------------------------------------------
const unsigned long long tab_num = 6;
const unsigned long long tab_front_module[tab_num] = {0, 3, 3, 2, 1, 0};
const unsigned long long tab_front_channel[tab_num] = {16, 0, 16, 16, 0, 0};
const unsigned long long tab_front_time_channel[tab_num] = {0, 80, 64, 48, 32, 16};
const unsigned long long tab_back_module[tab_num] = {1, 4, 4, 2, 2, 1};
const unsigned long long tab_back_channel[tab_num] = {24, 8, 0, 0, 8, 16};

const unsigned long long vtaf_num = 2;
const unsigned long long vtaf_front_module[vtaf_num] = {0, 0};
const unsigned long long vtaf_front_channel[vtaf_num] = {16, 0};
const unsigned long long vtaf_back_module[vtaf_num] = {1, 1};
const unsigned long long vtaf_back_channel[vtaf_num] = {8, 0};

//-----------------------------------------------------------------------------
//								detector config
//-----------------------------------------------------------------------------
const unsigned long long ppac_num = 3;

const std::vector<double> t0_thickness{
	1010.0, 1504.0, 1501.0, 1534.0, 1532.0, 1540.0
};

}

#endif 		// __GLOBAL_DEF_H__