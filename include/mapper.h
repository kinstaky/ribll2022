#ifndef __MAPPER_H__
#define __MAPPER_H__

#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/dssd_event.h"
#include "include/event/ppac_event.h"
#include "include/event/tof_event.h"

namespace ribll {



class XiaMapper {
public:

	/// @brief constructor
	///
	/// @param[in] run run number
	///
	XiaMapper(unsigned int run);


	/// @brief destructor
	///
	virtual ~XiaMapper();


	/// @brief abstract mapping function
	///
	/// @returns 0 if success, -1 on failure
	///
	virtual int Map() = 0;


	/// @brief initailize the input tree
	///
	/// @param[in] file_name name of input file
	/// @returns pointer to input tree if success, or nullptr otherwise
	///
	TTree* Initialize(const char *file_name);


	/// @brief calculate timestamp with raw timestamp and sampling rate
	///
	/// @param[in] rate sampling rate, 100, 250 or 500
	/// @param[in] ts timestamp read from file
	/// @returns real timestamp in nanosecond, or 0 for invalid rate
	///
	Long64_t CalculateTimestamp(Short_t rate, Long64_t ts);


	/// @brief calculate real time in nanosecond
	///
	/// @param[in] rate sampling rate, 100, 250 or 500
	/// @param[in] timestamp timestamp in nanosecond
	/// @param[in] cfd cfd value read from file
	/// @param[in] cfds cfd source value read from file
	/// @param[in] cfdft cfd forced trigger bit read from file
	/// @returns time value with cfd correction in nanosecond, 0 for invalid rate
	///
	Double_t CalculateTime(
		Short_t rate,
		Long64_t timestamp,
		Short_t cfd,
		Short_t cfds,
		Bool_t cfdft
	);


	/// @brief create residual tree
	///
	/// @param[in] crate name
	/// @returns index of tree
	///
	size_t CreateResidualTree(const char *name);


	/// @brief create output tree for detectors
	///
	/// @param[in] name detector name
	/// @returns index of tree
	///
	size_t CreateOutputTree(const char *name);


	/// @brief create output tree for triggers
	///
	/// @param[in] name detector name
	/// @returns index of tree
	///
	size_t CreateTriggerTree(const char *name);


	/// @brief fill event into output tree by index
	///
	/// @param[in] index index of tree, get from CreateOutputTree
	///
	inline void FillTree(size_t index) {
		opfs_[index]->cd();
		opts_[index]->Fill();
	}

protected:
	unsigned int run_;

	// input files
	std::vector<TFile*> ipfs_;

	// input data from decode root file
	Short_t rate_;
	Short_t sid_;
	Short_t ch_;
	Long64_t ts_;
	Short_t cfd_;
	Short_t cfds_;
	Bool_t cfdft_;
	UShort_t raw_energy_;

	// output files
	std::vector<TFile*> opfs_;
	// output trees
	std::vector<TTree*> opts_;

	// output data for output tree
	UShort_t detector_index_;
	UShort_t side_;
	UShort_t strip_;
	Long64_t timestamp_;
	Double_t time_;
	Double_t energy_;
};


class Crate0Mapper : public XiaMapper {
public:

	/// @brief constructor
	///
	/// @param[in] run run number
	///
	Crate0Mapper(unsigned int run);


	/// @brief default destructor
	///
	virtual ~Crate0Mapper() = default;


	/// @brief mapping function
	///
	/// @returns 0 for success, -1 otherwise
	///
	virtual int Map();
};



class Crate1Mapper : public XiaMapper {
public:

	/// @brief constructor
	///
	/// @param[in] run run number
	///
	Crate1Mapper(unsigned int run);


	/// @brief default destructor
	///
	virtual ~Crate1Mapper() = default;


	/// @brief mapping function
	///
	/// @returns 0 for success, -1 otherwise
	///
	virtual int Map();
};

class Crate2Mapper : public XiaMapper {
public:

	/// @brief constructor
	///
	/// @param[in] run run number
	///
	Crate2Mapper(unsigned int run);


	/// @brief default destructor
	///
	virtual ~Crate2Mapper() = default;


	/// @brief mapping function
	///
	/// @returns 0 for success, -1 otherwise
	///
	virtual int Map();
};


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
	/// @returns 0 for success, -1 otherwise
	///
	int Map();

private:
	// run number
	int run_;

	// input data
	Int_t adc_[5][32];
	Int_t madc_[2][32];
	Int_t gdc_[2][128][5];
	Int_t gmulti_[2][128];

	// output files
	std::vector<TFile*> opfs_;
	// output trees
	std::vector<TTree*> opts_;

	// output data
	PpacFundamentalEvent ppac_event_;
	DssdFundamentalEvent dssd_event_;
	TofFundamentalEvent tof_event_;
	long long align_time_;



	/// @brief initailize the input tree
	/// @param[in] file_name name of input file
	/// @returns pointer to input tree if success, or nullptr otherwise
	///
	TTree* Initialize(const char *file_name);


	/// @brief create ppac output tree
	/// @returns index of tree
	///
	size_t CreatePPACTree();


	/// @brief create ADSSD output tree
	/// @param[in] name detector name
	/// @returns index of tree
	///
	size_t CreateADSSDTree(const char *name);



	/// @brief create ToF output tree
	/// @returns index of tree
	///
	size_t CreateTofTree();


	/// @brief fill data into output tree by index
	/// @param[in] index index of tree to fill
	///
	inline void FillTree(size_t index) {
		opfs_[index]->cd();
		opts_[index]->Fill();
	}
};


}

#endif 			// __MAPPER_H__