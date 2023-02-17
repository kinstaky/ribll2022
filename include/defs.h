#ifndef __GLOBAL_DEF_H__
#define __GLOBAL_DEF_H__

namespace ribll {

//-----------------------------------------------------------------------------
//									data path
//-----------------------------------------------------------------------------
// decode data
const char* const kCrate1Path = "/data/d1/RIBLL_2022_C/DecodeFile/Crate1/";
const char* const kCrate1FileName = "c1data";

const char* const kCrate2Path = "/data/d1/RIBLL_2022_C/DecodeFile/Crate2/";
const char* const kCrate2FileName = "c2data";

const char* const kCrate3Path = "/data/d1/RIBLL_2022_C/DecodeFile/Crate3/";
const char* const kCrate3FileName = "c3data";

const char* const kCrate4Path = "/data/d1/RIBLL_2022_C/DecodeFile/VME/";
const char* const kCrate4FileName = "vmedata";

// generate data
const char* const kGenerateDataPath = "/data/d1/pwl/ribll2022/";
const char* const kMappingDir = "mapping/";
const char* const kDecodeDir = "decode/";
const char* const kAlignDir = "align/";
const char* const kFundamentalDir = "fundamental/";
const char* const kSingleSideDir = "single-side/";
const char* const kCorrelationDir = "correlation/";
const char* const kNormalizeDir = "normalize/";
const char* const kMergedDir = "merge/";
const char* const kTelescopeDir = "telescope/";


//-----------------------------------------------------------------------------
//								xia vme alignment
//-----------------------------------------------------------------------------
const int kScalerIndex = 4;
const int kScalerPeriod = 200;

//-----------------------------------------------------------------------------
//								vme module config
//-----------------------------------------------------------------------------
const unsigned long long tab_num = 6;
const unsigned long long tab_front_module[tab_num] = {0, 0, 1, 2, 3, 3};
const unsigned long long tab_front_channel[tab_num] = {16, 0, 0, 16, 16, 0};
const unsigned long long tab_back_module[tab_num] = {1, 1, 2, 2, 4, 4};
const unsigned long long tab_back_channel[tab_num] = {24, 16, 8, 0, 0, 8};

const unsigned long long vtaf_num = 2;
const unsigned long long vtaf_front_module[vtaf_num] = {0, 0};
const unsigned long long vtaf_front_channel[vtaf_num] = {16, 0};
const unsigned long long vtaf_back_module[vtaf_num] = {1, 1};
const unsigned long long vtaf_back_channel[vtaf_num] = {8, 0};

//-----------------------------------------------------------------------------
//								detector config
//-----------------------------------------------------------------------------
const unsigned long long ppac_num = 3;


}

#endif 		// __GLOBAL_DEF_H__