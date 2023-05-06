#ifndef __XIA_MAPPER_H__
#define __XIA_MAPPER_H__

#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "include/statistics/map_statistics.h"

namespace ribll {

class XiaMapper {
public:

	/// @brief constructor
	/// @param[in] run run number
	///
	XiaMapper(unsigned int run);


	/// @brief destructor
	///
	virtual ~XiaMapper();


	/// @brief abstract mapping function
	/// @param[in] threshold whether to use threshold
	/// @returns 0 if success, -1 on failure
	///
	virtual int Map(bool threshold = true) = 0;

protected:
	unsigned int run_;

	// input files
	std::vector<TFile*> ipfs_;

	// input data from decode root file
	short rate_;
	short sid_;
	short ch_;
	long long ts_;
	short cfd_;
	short cfds_;
	bool cfdft_;
	unsigned short raw_energy_;

	// output files
	std::vector<TFile*> opfs_;
	// output trees
	std::vector<TTree*> opts_;

	// output data for output tree
	unsigned short detector_index_;
	unsigned short side_;
	unsigned short strip_;
	long long timestamp_;
	double time_;
	double energy_;
	long long decode_entry_;


	/// @brief initailize the input tree
	///
	/// @param[in] file_name name of input file
	/// @returns pointer to input tree if success, or nullptr otherwise
	///
	TTree* Initialize(const char *file_name);


	/// @brief calculate timestamp with raw timestamp and sampling rate
	/// @param[in] rate sampling rate, 100, 250 or 500
	/// @param[in] ts timestamp read from file
	/// @returns real timestamp in nanosecond, or 0 for invalid rate
	///
	long long CalculateTimestamp(short rate, long long ts);


	/// @brief calculate real time in nanosecond
	/// @param[in] rate sampling rate, 100, 250 or 500
	/// @param[in] timestamp timestamp in nanosecond
	/// @param[in] cfd cfd value read from file
	/// @param[in] cfds cfd source value read from file
	/// @param[in] cfdft cfd forced trigger bit read from file
	/// @returns time value with cfd correction in nanosecond,
	/// 	0 for invalid rate
	///
	double CalculateTime(
		short rate,
		long long timestamp,
		short cfd,
		short cfds,
		bool cfdft
	);


	/// @brief create residual tree
	/// @param[in] crate name
	/// @returns index of tree
	///
	size_t CreateResidualTree(const char *name);


	/// @brief create output tree for detectors
	/// @param[in] name detector name
	/// @param[in] threshold whether to add nc- to file name (nc == no cut)
	/// @returns index of tree
	///
	size_t CreateOutputTree(const char *name, bool threshold);


	/// @brief create output tree for triggers
	/// @param[in] name detector name
	/// @param[in] threshold whether to add nc- to file name (nc == no cut)
	/// @returns index of tree
	///
	size_t CreateTriggerTree(const char *name, bool threshold);


	/// @brief fill event into output tree by index
	/// @param[in] index index of tree, get from CreateOutputTree
	///
	inline void FillTree(size_t index) {
		opfs_[index]->cd();
		opts_[index]->Fill();
	}
};

}		// namespace ribll

#endif 			// __XIA_MAPPER_H__