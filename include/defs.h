#ifndef __GLOBAL_DEF_H__
#define __GLOBAL_DEF_H__

#include <vector>
#include <cmath>

namespace ribll {

//-----------------------------------------------------------------------------
//									data path
//-----------------------------------------------------------------------------
// decode data
const char* const kCrate0Path = "/data/d1/RIBLL_2022_C/DecodeFile/Crate1/";
// const char* const kCrate0Path = "/data/d1/pwl/ribll2022/decode/";
const char* const kCrate0FileName = "c1data";

const char* const kCrate1Path = "/data/d1/RIBLL_2022_C/DecodeFile/Crate2/";
const char* const kCrate1FileName = "c2data";

const char* const kCrate2Path = "/data/d1/RIBLL_2022_C/DecodeFile/Crate3/";
const char* const kCrate2FileName = "c3data";

const char* const kCrate3Path = "/data/d1/RIBLL_2022_C/DecodeFile/VME/";
const char* const kCrate3FileName = "vmedata";

// generate data
const char* const kGenerateDataPath = "/data/d1/pwl/ribll2022/";
const char* const kMappingDir = "map/";
const char* const kDecodeDir = "decode/";
const char* const kAlignDir = "align/";
const char* const kFundamentalDir = "fundamental/";
const char* const kTraceDir = "trace/";
const char* const kTimeDir = "time/";
const char* const kBeamDir = "beam/";
const char* const kNormalizeDir = "normalize/";
const char* const kMergeDir = "merge/";
const char* const kTelescopeDir = "telescope/";
const char* const kCalibrationDir = "calibration/";
const char* const kParticleIdentifyDir = "pid/";
const char* const kParticleDir = "particle/";
const char* const kChannelDir = "channel/";
const char* const kSpectrumDir = "spectrum/";

const char* const kCheckDir = "check/";
const char* const kShowDir = "show/";
const char* const kEnergyCalculateDir = "energy_calculate/";
const char* const kFilterDir = "filter/";
const char* const kSummaryDir = "summary/";
const char* const kSimulateDir = "simulate/";
const char* const kInformationDir = "info/";
const char* const kOptimizeDir = "opt/";
const char* const kHoleDir = "hole/";

//-----------------------------------------------------------------------------
//								xia vme alignment
//-----------------------------------------------------------------------------
const int kScalerIndex = 4;
const int kScalerPeriod = 200;

//-----------------------------------------------------------------------------
//								vme module config
//-----------------------------------------------------------------------------
const unsigned int tab_num = 6;
const unsigned int tab_front_module[tab_num] = {0, 3, 3, 2, 1, 0};
const unsigned int tab_front_channel[tab_num] = {16, 0, 16, 16, 0, 0};
const unsigned int tab_front_time_channel[tab_num] = {0, 80, 64, 48, 32, 16};
const unsigned int tab_back_module[tab_num] = {1, 4, 4, 2, 2, 1};
const unsigned int tab_back_channel[tab_num] = {24, 8, 0, 0, 8, 16};

const unsigned int vtaf_num = 2;
const unsigned int vtaf_front_module[vtaf_num] = {0, 0};
const unsigned int vtaf_front_channel[vtaf_num] = {16, 0};
const unsigned int vtaf_back_module[vtaf_num] = {1, 1};
const unsigned int vtaf_back_channel[vtaf_num] = {8, 0};

//-----------------------------------------------------------------------------
//								detector config
//-----------------------------------------------------------------------------

const std::vector<double> t0_thickness {
	1010.0, 1504.0, 1501.0, 1534.0, 1532.0, 1540.0
};
constexpr double t0z[3] = {100.0, 111.76, 123.52};

constexpr unsigned int ppac_num = 3;
constexpr double ppac_xz[ppac_num] = {-695.2, -454.2, -275.2};
constexpr double ppac_yz[ppac_num] = {-689.2, -448.2, -269.2};
constexpr double all_ppac_xz[4] = {-695.2, -633.7, -454.2, -275.2};
constexpr double all_ppac_yz[4] = {-689.2, -627.7, -448.2, -269.2};
constexpr int ppac_change_run = 717;

// thickness of tafd
constexpr double tafd_thickness[6] = {166.0, 154.0, 158.0, 162.0, 150.0, 164.0};


//-----------------------------------------------------------------------------
//								calibration parameters
//-----------------------------------------------------------------------------

constexpr double t0_param[6][2] = {
	{0.0553516, 0.00532019},
	{-0.12591, 0.00632308},
	{0.552785, 0.00579009},
	{0.837779, 0.00233567},
	{-0.306592, 0.00221028},
	{3.0818, 0.00235991}
};

constexpr double csi_param[12][3] = {
	{216.579, 0.97, -23.4389},
	{214.067, 0.96, -26.1516},
	{200.308, 1.02, -155.26},
	{276.858, 0.96, -468.231},
	{233.155, 1.02, -211.293},
	{327.988, 0.96, -988.938},
	{318.893, 0.96, -758.227},
	{269.045, 0.96, -222.745},
	{291.919, 0.96, -346.013},
	{241.453, 0.96, 100.75},
	{268.52, 0.96, -196.279},
	{234.844, 1.02, -188.712}
};

constexpr double power_csi_param[12][3] = {
	{116.414, 1.12316, 118.003},
	{115.785, 1.11126, 87.6736},
	{123.988, 1.1368, 330.016},
	{183.623, 1.05357, 47.4882},
	{202.189, 1.05736, -21.0135},
	{190.769, 1.06031, 53.1218},
	{187.891, 1.07026, -27.4448},
	{167.966, 1.08199, 28.1004},
	{214.845, 1.06909, -13.6187},
	{235.779, 1.0131, -211.072},
	{156.749, 1.09095, 29.8454},
	{183.084, 1.08151, -84.3011}
};

//-----------------------------------------------------------------------------
//								physical constants
//-----------------------------------------------------------------------------
constexpr double pi = 3.14159265359;

constexpr double u = 931.494;
constexpr double mass_n = u * 1.00866491588;
constexpr double mass_1h = u * 1.0072764520;
constexpr double mass_2h = u * 2.0135531980;
constexpr double mass_4he = u * 4.0015060943;
constexpr double mass_9be = u * 9.0099887420;
constexpr double mass_10be = u * 10.0113403769;
constexpr double mass_14c = u * 13.9999505089;
constexpr double mass_15c = u * 15.0073077289;


/// @brief get momentum from kinematic energy
/// @param[in] mass mass of particle
/// @param[in] kinematic kinematic energy of particle
/// @returns momentum of particle
///
inline double MomentumFromKinetic(double mass, double kinematic) {
	return sqrt((2.0 * mass + kinematic) * kinematic);
}


}	// ribll

#endif 		// __GLOBAL_DEF_H__