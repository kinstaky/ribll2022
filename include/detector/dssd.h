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
	//							normalize
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


	/// @brief generate filter flag for iteration mode
	/// @param[in] iteration iteration number
	/// @returns 0 if success, -1 otherwise
	///
	virtual int NormalizeFilter(int iteration);

	//-------------------------------------------------------------------------
	//							merge
	//-------------------------------------------------------------------------

	/// @brief merge adjacent event in the same side and merge events of two sides
	/// @param[in] energy_diff tolerant energy relateive difference
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Merge(double energy_diff);

protected:

	//-------------------------------------------------------------------------
	//						geometry function
	//-------------------------------------------------------------------------

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


	/// @brief normalize one side
	/// @param[in] chain input TChain
	/// @param[in] side side to normalize, 0 for front, 1 for back
	/// @param[in] ref_strip reference strip from the other side
	/// @param[in] iteration iteration mode
	///
	virtual int SideNormalize(
		TChain *chain,
		size_t side,
		size_t ref_strip,
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
		return norm_params_[side][strip][0]
			+ norm_params_[side][strip][1] * energy;
	}


	/// @brief read cut from file
	/// @param[in] name file name
	///
	std::unique_ptr<TCutG> ReadCut(const std::string &name) const;

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
	double norm_params_[2][64][2];
};


inline double RelativeDifference(double x, double y) {
	return fabs((x - y) / (x + y));
}


}		// namespace ribll

#endif		// __DSSD_H__