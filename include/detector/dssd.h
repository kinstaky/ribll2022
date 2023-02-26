#ifndef __DSSD_H__
#define __DSSD_H__

#include <string>

#include "include/detector/detector.h"

namespace ribll {

class Dssd : public Detector {
public:
	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	///
	Dssd(unsigned int run, const std::string &name);


	/// @brief default destructor
	///
	virtual ~Dssd() = default;


	/// @brief statistics of matching trigger
	///
	class MatchTriggerStatistics {
	public:
		long long total_events;
		long long match_events;
		long long oversize_events;

		/// @brief constructor
		/// @param[in] total total events
		/// 
		MatchTriggerStatistics(long long total);


		/// @brief overloaded operator<< function, output statistics
		/// @param[in] os ostream
		/// @param[in] statistics output object
		/// @returns ostream
		///
		friend std::ostream& operator<<(
			std::ostream &os,
			const MatchTriggerStatistics &statisics
		);
	};


	/// @brief match xia main trigger and build events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(double window_left, double window_right);


	// /// @brief get front strip number
	// ///
	// /// @returns front strip number
	// ///
	// virtual inline size_t FrontStrip() override
	// {
	// 	return 32;
	// }


	// /// @brief get back strip number
	// ///
	// /// @returns back strip number
	// ///
	// virtual inline size_t BackStrip() override
	// {
	// 	return 32;
	// }


	// /// @brief match xia main trigger and build events
	// /// @param[in] window_left left edge of match window
	// /// @param[in] window_right right edge of match window
	// /// @returns 0 if success, -1 otherwise
	// ///
	// virtual int MatchTrigger(double window_left, double window_right);

// private:
// 	/// @brief read tirgger from root file
// 	/// @param[out] trigger_times list of trigger time
// 	/// @returns 0 if success, -1 otherwise
// 	///
// 	int ReadTriggerTimes(std::vector<double> &trigger_times);


// 	/// @brief setup output tree of correlated events
// 	///
// 	/// @param[in] tree output correlated tree to setup
// 	///
// 	virtual void SetupOutputCorrelatedTree(TTree *tree) override;


// 	/// @brief read correlated events from root file
// 	///
// 	/// @returns 0 for success, -1 otherwise
// 	///
// 	virtual void SetupInputCorrelatedTree(TTree *tree) override;
};


class T0d1 : public Dssd {
public:

	/// @brief constructor
	///
	/// @param[in] run run number
	///
	T0d1(unsigned int run);


	/// @brief default destructor
	///
	virtual ~T0d1() = default;


	// /// @brief get front strip number
	// ///
	// /// @returns front strip number
	// ///
	// inline size_t FrontStrip() override {
	// 	return 64;
	// }


	// /// @brief get back strip number
	// ///
	// /// @returns back strip number
	// ///
	// virtual inline size_t BackStrip() override {
	// 	return 64;
	// }

// private:
// 	/// @brief check whether the energy of correlated event can be use in normalization
// 	/// for back strip
// 	///
// 	/// @param[in] correlation corrlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns true if pass, or false otherwise
// 	///
// 	virtual bool NormalizeFrontEnergyCheck(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) override;


// 	/// @brief check whether the energy of correlated event can be use in normalization
// 	/// for front strip
// 	///
// 	/// @param[in] correlation correlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns ture if pass, or false otherwise
// 	///
// 	virtual bool NormalizeBackEnergyCheck(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) override;
};


class T0d2 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	///
	T0d2(unsigned int run);


	/// @brief default destructor
	///
	virtual ~T0d2() = default;


// 	/// @brief get front strip number
// 	///
// 	/// @returns front strip number
// 	///
// 	inline size_t FrontStrip() override
// 	{
// 		return 32;
// 	}


// 	/// @brief get back strip number
// 	///
// 	/// @returns back strip number
// 	///
// 	inline size_t BackStrip() override
// 	{
// 		return 32;
// 	}


// private:
// 	/// @brief check whether the energy of correlated event can be use in normalization
// 	/// for back strip
// 	///
// 	/// @param[in] correlation corrlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns true if pass, or false otherwise
// 	///
// 	virtual bool NormalizeFrontEnergyCheck(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) override;


// 	/// @brief check whether the energy of correlated event can be use in normalization
// 	/// for front strip
// 	///
// 	/// @param[in] correlation correlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns ture if pass, or false otherwise
// 	///
// 	virtual bool NormalizeBackEnergyCheck(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) override;
};


class T0d3 : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	///
	T0d3(unsigned int run);


	/// @brief default destructor
	///
	virtual ~T0d3() = default;


// 	/// @brief get front strip number
// 	///
// 	/// @returns front strip number
// 	///
// 	inline size_t FrontStrip() override {
// 		return 32;
// 	}


// 	/// @brief get back strip number
// 	///
// 	/// @returns back strip number
// 	///
// 	virtual inline size_t BackStrip() override {
// 		return 32;
// 	}

// private:
// 	/// @brief check whether the energy of correlated event can be use in normalization
// 	/// for back strip
// 	///
// 	/// @param[in] correlation corrlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns true if pass, or false otherwise
// 	///
// 	virtual bool NormalizeFrontEnergyCheck(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) override;


// 	/// @brief check whether the energy of correlated event can be use in normalization
// 	/// for front strip
// 	///
// 	/// @param[in] correlation correlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns ture if pass, or false otherwise
// 	///
// 	virtual bool NormalizeBackEnergyCheck(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) override;
};


}


#endif				// __DSSD_H__