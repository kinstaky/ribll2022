#ifndef __CRATE3_MAPPER_H__
#define __CRATE3_MAPPER_H__

#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/dssd_event.h"
#include "include/event/ppac_event.h"
#include "include/event/tof_event.h"

namespace ribll {

class Crate3Mapper {
public:
	/// @brief constructor
	/// @param[in] run run number
	///
	Crate3Mapper(unsigned int run);


	/// @brief default destructor
	///
	virtual ~Crate3Mapper() = default;


	/// @brief mapping function
	/// @param[in] independent independent to XIA
	/// @returns 0 for success, -1 otherwise
	///
	int Map(bool independent = false);

private:

	/// @brief initailize the input tree
	/// @param[in] file_name name of input file
	/// @returns pointer to input tree if success, or nullptr otherwise
	///
	TTree* Initialize(const char *file_name);


	/// @brief create ppac output tree
	/// @param[in] independent independent to XIA
	/// @returns index of tree
	///
	size_t CreatePPACTree(bool independent = false);


	/// @brief create ADSSD output tree
	/// @param[in] name detector name
	/// @param[in] independent independent to XIA
	/// @returns index of tree
	///
	size_t CreateADSSDTree(const char *name, bool independent = false);



	/// @brief create ToF output tree
	/// @param[in] independent independent to XIA
	/// @returns index of tree
	///
	size_t CreateTofTree(bool independent = false);


	/// @brief fill data into output tree by index
	/// @param[in] index index of tree to fill
	///
	inline void FillTree(size_t index) {
		opfs_[index]->cd();
		opts_[index]->Fill();
	}


	// run number
	int run_;

	// input data
	int adc_[5][32];
	int madc_[2][32];
	int gdc_[2][128][5];
	int gmulti_[2][128];

	// output files
	std::vector<TFile*> opfs_;
	// output trees
	std::vector<TTree*> opts_;

	// output data
	PpacFundamentalEvent ppac_event_;
	DssdFundamentalEvent dssd_event_;
	TofFundamentalEvent tof_event_;
	long long align_time_;
	int align_gmulti_[128];
	int align_gdc_[128][5];
};

}		// namespace ribll

#endif 		// __CRATE3_MAPPER_H__