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


//-----------------------------------------------------------------------------
//								xia vme alignment
//-----------------------------------------------------------------------------
const int kScalerIndex = 4;
const int kScalerPeriod = 200;

}

#endif 		// __GLOBAL_DEF_H__