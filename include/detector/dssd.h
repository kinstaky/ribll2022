#ifndef __DSSD_H__
#define __DSSD_H__

#include <string>

#include <TChain.h>
#include <TCutG.h>
#include <TMath.h>
#include <Math/Vector3D.h>

#include "include/detector/detector.h"
#include "include/event/dssd_event.h"
#include "include/statistics/merge_statistics.h"
#include "include/statistics/merge_case_statistics.h"

namespace ribll {

struct NormalizeInfo {
	int side;
	unsigned short ref_start;
	unsigned short ref_end;
	unsigned short norm_start;
	unsigned short norm_end;
};

class Dssd : public Detector {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	/// @param[in] tag trigger tag
	///
	Dssd(
		unsigned int run,
		const std::string &name,
		const std::string &tag
	);


	/// @brief default destructor
	///
	virtual ~Dssd() = default;


	//-------------------------------------------------------------------------
	//								gemometry
	//-------------------------------------------------------------------------

	/// @brief get front strip number
	/// @returns front strip number
	///
	virtual inline size_t FrontStrip() const {
		return 32;
	}


	/// @brief get back strip number
	/// @returns back strip number
	///
	virtual inline size_t BackStrip() const {
		return 32;
	}


	/// @brief return strip number based on side
	/// @param[in] side 0 for front side, 1 for back side
	/// @returns strip number
	///
	virtual inline size_t Strip(int side) const {
		return side == 0 ? FrontStrip() : BackStrip();
	}


	/// @brief calculate the position from strip index
	/// @param[in] front_strip front strip
	/// @param[in] back_strip front strip
	/// @returns vector point to the position
	///
	virtual ROOT::Math::XYZVector CalculatePosition(
		double front_strip,
		double back_strip
	) const;

	//-------------------------------------------------------------------------
	//							match trigger
	//-------------------------------------------------------------------------

	/// @brief match xia main trigger and build events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		double window_left,
		double window_right
	) override;


	/// @brief extract trigger with detector events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ExtractTrigger(
		double window_left,
		double window_right
	) override;


	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	/// @brief normalize dssd
	/// @param[in] end_run end of run to chain, inclusive
	/// @param[in] iteration iteration mode, default is 0 for normal mode
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Normalize(
		unsigned int end_run,
		int iteration = 0
	);


	/// @brief show normalize result
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ShowNormalize();


	/// @brief wirte normalized energy to root file
	/// @param[in] iteration iteration mode
	/// @returns 0 if success, -1 otherwise
	///
	int NormalizeResult(int iteration = 0);

	//-------------------------------------------------------------------------
	//								merge
	//-------------------------------------------------------------------------


	/// @brief cut beam or events under threshold
	/// @returns 0 if success, -1 otherwise
	///
	virtual int CutBeamThreshold();

	/// @brief merge adjacent event in the same side and merge events of two sides
	/// @param[in] energy_diff tolerant energy relateive difference
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Merge(double energy_diff);

protected:

	/// @brief Fill merge event from normalize result event, version 2
	/// @param[in] event normalize result event
	/// @param[in] energy_diff energy difference tolerance
	/// @param[out] merge merge event to fill
	/// @returns merge hit
	int FillMergeEvent2(
		const DssdFundamentalEvent &event,
		double energy_diff,
		DssdMergeEvent &merge
	);


	/// @brief search front ajacent strip to specific strip
	/// @param[in] event normalize result event
	/// @param[in] strip strip to search
	/// @param[in] flag used strip flag
	/// @returns index of adjacent strip if found, -1 otherwise
	///
	int SearchFrontAdjacentStrips(
		const DssdFundamentalEvent &event,
		unsigned short strip,
		unsigned short flag
	);


	/// @brief search back ajacent strip to specific strip
	/// @param[in] event normalize result event
	/// @param[in] strip strip to search
	/// @param[in] flag used strip flag
	/// @returns index of adjacent strip if found, -1 otherwise
	///
	int SearchBackAdjacentStrips(
		const DssdFundamentalEvent &event,
		unsigned short strip,
		unsigned short flag
	);


	/// @brief sort particles in merge event by energy
	/// @param[inout] merge merge event to sort
	///
	void SortMergeEvent(DssdMergeEvent &merge);


	//-------------------------------------------------------------------------
	//								time
	//-------------------------------------------------------------------------

public:
	/// @brief analyze time
	/// @returns 0 if success, -1 otherwise
	///
	virtual int AnalyzeTime();


	/// @brief normalize time
	/// @returns 0 if success, -1 otherwise
	///
	virtual int NormalizeTime();


	/// @brief filter curve in time-energy histogram
	/// @returns 0 if success, -1 otherwise
	///
	virtual int FilterTimeCurve();


	/// @brief fit time curve in time-energy histogram
	/// @returns 0 if success, -1 otherwise
	///
	virtual int FitTimeCurve();


	/// @brief check whether time curve is appropriate
	/// @param[in] condition 0-single hit, 1-double hit small, 2-double hit big
	/// @param[in] side front (0) or back (1)
	/// @param[in] strip strip number
	/// @param[in] energy normalized energy
	/// @param[in] time normalized time
	/// @returns true if pass time check, false otherwise
	virtual bool CheckTime(
		int condition,
		size_t side,
		unsigned short strip,
		double energy,
		double time
	);

protected:

	//-------------------------------------------------------------------------
	//						normalize function
	//-------------------------------------------------------------------------

	/// @brief read normalize parameters form file
	/// @param iteration itertation mode,
	/// 	-1 to read the file without tag
	/// @returns 0 if success, -1 otherwise
	///
	int ReadNormalizeParameters(int iteration = -1);


	/// @brief write normalize parameters to file
	/// @param iteration itertation mode,
	///		-1 for overwriting parameter file without tag
	/// @returns 0 if success, -1 otherwise
	///
	int WriteNormalizeParameters(int iteration = -1);


	/// @brief normalize some strips in reference some strips in the other side
	/// @param[in] chain input TChain
	/// @param[in] side side to normalize, 0 for front, 1 for back
	/// @param[in] ref_start start strip to reference, inclusive
	/// @param[in] ref_end end strip to reference, exclusive
	/// @param[in] norm_start start strip to normalize, inclusive
	/// @param[in] norm_end end strip to normalize, exclusive
	/// @param[in] iteration iteration mode
	/// @returns 0 if success, -1 otherwise
	///
	virtual int StripsNormalize(
		TChain *chain,
		size_t side,
		size_t ref_start,
		size_t ref_end,
		size_t norm_start,
		size_t norm_end,
		int iteration
	);


	/// @brief normalize both sides, the true normalize
	/// @param[in] chain TChain of input events
	/// @param[in] iteration iteration mode
	/// @returns 0 if success, -1 otherwise
	///
	virtual int NormalizeSides(
		TChain *chain,
		int iteration
	);


	/// @brief get normalized energy based on parameters
	/// @param[in] side 0 for front, 1 for back
	/// @param[in] strip strip index
	/// @param[in] energy energy
	/// @returns normaized energy
	///
	inline double NormEnergy(int side, int strip, double energy) {
		return norm_params_[side][strip][0] / (energy + 1e-5)
			+ norm_params_[side][strip][1]
			+ norm_params_[side][strip][2] * energy
			+ norm_params_[side][strip][3] * energy * energy;
	}


	/// @brief read cut from file
	/// @param[in] dir directory name
	/// @param[in] name file name
	/// @returns pointer to cut object
	///
	std::unique_ptr<TCutG> ReadCut(const std::string &dir, const std::string &name) const;


	//-------------------------------------------------------------------------
	//							time function
	//-------------------------------------------------------------------------

	/// @brief read time normalized parameters from file
	/// @returns 0 if success, -1 otherwise
	///
	int ReadNormalizeTimeParameters();


	/// @brief write time normalized parameters to file
	/// @returns 0 if success, -1 otherwise
	///
	int WriteNormalizeTimeParameters();


	/// @brief read time cuts
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ReadTimeCuts();


	//-------------------------------------------------------------------------
	//						geometry member
	//-------------------------------------------------------------------------

	ROOT::Math::XYZVector center_;
	std::pair<double, double> x_range_;
	std::pair<double, double> y_range_;

	//-------------------------------------------------------------------------
	//						normalize member
	//-------------------------------------------------------------------------

	// normalize parameters, first index is side,
	// second index is strip, third index is p0 and p1
	double norm_params_[2][64][4];
	// time normalize parameters, first index is side,
	// second index is strip
	double norm_time_params_[2][64];

	//-------------------------------------------------------------------------
	//							time memeber
	//-------------------------------------------------------------------------

	std::vector<std::unique_ptr<TCutG>> time_cuts_[3];

private:
	bool has_normalized_[2][64];
};


inline double RelativeDifference(double x, double y) {
	return fabs((x - y) / (x + y));
}


}		// namespace ribll

#endif		// __DSSD_H__