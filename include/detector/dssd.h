#ifndef __DSSD_H__
#define __DSSD_H__

#include <string>

#include <TChain.h>

#include "include/detector/detector.h"
#include "include/event/dssd_event.h"

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
	/// @param[in] length number of files to use
	/// @param[in] iteration in iteration mode?
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Normalize(
		unsigned int length,
		bool iteration
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


	/// @brief show normalize result
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ShowNormalize();


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
	//							normalize
	//-------------------------------------------------------------------------

	/// @brief read normalize parameters form file
	/// @returns 0 if success, -1 otherwise
	///
	int ReadNormalizeParameters();


	/// @brief write normalize parameters to file
	/// @returns 0 if success, -1 otherwise
	///
	int WriteNormalizeParameters();


	/// @brief normalize one side
	/// @param[in] chain input TChain
	/// @param[in] side side to normalize, 0 for front, 1 for back
	/// @param[in] ref_strip reference strip from the other side
	/// @param[in] iteration iteration mode
	///
	int SideNormalize(
		TChain *chain,
		size_t side,
		size_t ref_strip,
		bool iteration
	);


	/// @brief normalize both sides, the true normalize
	/// @param[in] chain TChain of input events
	/// @param[in] iteration iteration mode?
	/// @returns 0 if success, -1 otherwise
	///
	virtual int NormalizeSides(TChain *chain, bool iteration);


	/// @brief wirte normalized energy to root file
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	/// @returns 0 if success, -1 otherwise
	///
	int WriteNormalizeFiles(unsigned int run, const std::string &tag);


	/// @brief check whether energy is suitable for fitting
	/// @param[in] side side to normalize
	/// @param[in] event fundamental event
	/// @returns true if pass check, false not pass
	///
	virtual bool NormEnergyCheck(
		size_t side,
		const DssdFundamentalEvent &event
	) const;


	// normalize parameters, first index is side,
	// second index is strip, third index is p0 and p1
	double norm_params_[2][64][2];
};


inline double RelativeDifference(double x, double y) {
	return fabs((x - y) / (x + y));
}


}		// namespace ribll

#endif		// __DSSD_H__