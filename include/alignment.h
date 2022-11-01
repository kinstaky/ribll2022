#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include <vector>
#include <map>

#include <TGraph.h>

namespace ribll {

class Alignment {
public:

	Alignment(
		int run,
		size_t group_num,
		double search_window,
		double search_low_bound,
		double search_high_bound
	);

	virtual ~Alignment();

	virtual int Align();

	inline void SetVerbose(bool verbose = true) {
		verbose_ = verbose;
	}

private:
	// run number
	int run_;

	// group number
	int group_num_;
	// search window size(step), low and high bound
	double search_window_;
	double search_low_bound_;
	double search_high_bound_;

	// xia data map
	std::vector<Long64_t> xia_times_;
	// vme time list
	std::vector<Long64_t> vme_times_;

	// print processing information
	bool verbose_;

	int ReadXiaTime();

	int ReadVmeTime();

	TGraph* GroupAlignment();

	int BuildResult(Double_t *calibration_param);
};

}


#endif			// __ALIGNMENT_H__