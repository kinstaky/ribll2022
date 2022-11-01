#ifndef __GLOBAL_DEF_H__
#define __GLOBAL_DEF_H__

namespace ribll {

//-----------------------------------------------------------------------------
//									data path
//-----------------------------------------------------------------------------
// decode data
const char *kCrate1Path = "/mnt/crate1/decode/";
const char *kCrate1FileName = "c1data";

const char *kCrate2Path = "/mnt/crate2/decode/";
const char *kCrate2FileName = "c2data";

const char *kCrate3Path = "/mnt/crate3/decode/";
const char *kCrate3FileName = "c3data";

const char *kCrate4Path = "/mnt/vme/rootfile/";
const char *kCrate4FileName = "vmedata";

// generate data
const char *kGenerateDataPath = "/home/test/Analysis/ribll2022/data/";
const char *kMappingDir = "mapping/";
const char *kDecodeDir = "decode/";
const char *kAlignDir = "align/";
const char *kCorrelationDir = "correlation/";


//-----------------------------------------------------------------------------
//								xia vme alignment
//-----------------------------------------------------------------------------
const int kScalerIndex = 4;
const int kScalerPeriod = 200;


//-----------------------------------------------------------------------------
//								vme module config
//-----------------------------------------------------------------------------
const size_t tab_num = 6;
const size_t tab_front_module[tab_num] = {0, 0, 1, 2, 3, 3};
const size_t tab_front_channel[tab_num] = {16, 0, 0, 16, 16, 0};
const size_t tab_back_module[tab_num] = {1, 1, 2, 2, 4, 4};
const size_t tab_back_channel[tab_num] = {31, 23, 15, 7, 15, 7};

const size_t vtaf_num = 2;
const size_t vtaf_front_module[vtaf_num] = {0, 0};
const size_t vtaf_front_channel[vtaf_num] = {16, 0};
const size_t vtaf_back_module[vtaf_num] = {1, 1};
const size_t vtaf_back_channel[vtaf_num] = {15, 7};

}

#endif 		// __GLOBAL_DEF_H__