#ifndef __T0D1_H__
#define __T0D1_H__

#include "include/detector/dssd.h"

namespace ribll {

class T0d1 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	///
	T0d1(unsigned int run, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~T0d1() = default;


	//-------------------------------------------------------------------------
	//								gemometry
	//-------------------------------------------------------------------------

	/// @brief get front strip number
	/// @returns front strip number
	///
	inline size_t FrontStrip() const override {
		return 64;
	}


	/// @brief get back strip number
	/// @returns back strip number
	///
	virtual inline size_t BackStrip() const override {
		return 64;
	}


	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	/// @brief normalize dssd
	/// @param[in] end_run end of run to chain, inclusive
	/// @param[in] iteration iteration mode, default is 0 for normal mode
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Normalize(unsigned int end_run, int iteration = 0);

	//-------------------------------------------------------------------------
	//								merge
	//-------------------------------------------------------------------------

	/// @brief cut beam or events under threshold
	/// @returns 0 if success, -1 otherwise
	///
	int CutBeamThreshold() override;

	//-------------------------------------------------------------------------
	//								time
	//-------------------------------------------------------------------------

	/// @brief analyze time
	/// @returns 0 if success, -1 otherwise
	///
	virtual int AnalyzeTime() override;


	/// @brief filter curve in time-energy histogram
	/// @returns 0 if success, -1 otherwise
	///
	virtual int FilterTimeCurve() override;


	/// @brief fit time curve in time-energy histogram
	/// @returns 0 if success, -1 otherwise
	///
	virtual int FitTimeCurve() override;


	/// @brief check whether time curve is appropriate
	/// @param[in] condition 0-single hit, 1-double hit small, 2-double hit big
	/// @param[in] side front (0) or back (1)
	/// @param[in] strip strip number
	/// @param[in] energy normalized energy
	/// @param[in] time normalized time
	/// @returns true if pass time check, false otherwise
	///
	virtual bool CheckTime(
		int condition,
		size_t side,
		unsigned short strip,
		double energy,
		double time
	) override;


protected:

	//-------------------------------------------------------------------------
	//								geometry
	//-------------------------------------------------------------------------

	/// @brief calculate the position from strip index
	/// @param[in] front_strip front strip
	/// @param[in] back_strip front strip
	/// @returns vector point to the position
	///
	virtual ROOT::Math::XYZVector CalculatePosition(
		double front_strip,
		double back_strip
	) const override;


	//-------------------------------------------------------------------------
	//								normalize
	//-------------------------------------------------------------------------

	/// @brief normalize both sides, the true normalize
	/// @param[in] chain TChain of input events
	/// @param[in] iteration iteration mode
	/// @returns 0 if success, -1 otherwise
	///
	virtual int NormalizeSides(
		TChain *chain,
		int iteration
	) override;


	/// @brief read time cuts
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ReadTimeCuts() override;
};



}		// namespace ribll

#endif		// __T0D1_H__