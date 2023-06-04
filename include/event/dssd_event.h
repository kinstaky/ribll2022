#ifndef __DSSD_EVENT_H__
#define __DSSD_EVENT_H__

#include "include/event/event.h"

namespace ribll {

class DssdMapEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	/// @param[in] prefix prefix of variables for friend tree
	///
	virtual void SetupInput(
		TTree *tree,
		const std::string &prefix = ""
	) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	unsigned short side;
	unsigned short strip;
	bool cfd_flag;
	double time;
	double energy;
	long long decode_entry;
};


class DssdFundamentalEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	/// @param[in] prefix prefix of variables for friend tree
	///
	virtual void SetupInput(
		TTree *tree,
		const std::string &prefix = ""
	) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	/// @brief show whether the event is valid or not
	/// @returns true if valid, false if invalid
	///
	virtual inline bool Valid() const {
		return front_hit != 0 && back_hit != 0;
	}


	/// @brief make the event invalid
	///
	virtual inline void Nullify() {
		front_hit = back_hit = 0;
	}


	/// @brief sort events by energy
	///
	virtual void Sort();


	/// @brief erase an event
	/// @param[in] side 0-front, 1-back
	/// @param[in] index index of event to be removed
	///
	virtual void Erase(size_t side, size_t index);

	unsigned short front_hit;
	unsigned short back_hit;
	unsigned short cfd_flag;
	unsigned short front_strip[8];
	unsigned short back_strip[8];
	double front_time[8];
	double back_time[8];
	double front_energy[8];
	double back_energy[8];
	long long front_decode_entry[8];
	long long back_decode_entry[8];

private:

	/// @brief sort one side by energy start from index
	/// @param[in] side 0-front, 1-back
	/// @param[in] index event before this index is sorted
	///
	void SortSide(size_t side, unsigned short index);


	/// @brief swap two event in one side
	/// @param[in] side 0-front, 1-back
	/// @param[in] i first index to swap
	/// @param[in] j second index to swap
	///
	void Swap(size_t side, unsigned short i, unsigned short j);
};


class DssdTimeEvent : public Event {
public:
	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	/// @param[in] prefix prefix of variables for friend tree
	///
	virtual void SetupInput(
		TTree *tree,
		const std::string &prefix = ""
	) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;

	unsigned short front_hit;
	unsigned short back_hit;
	int front_time_flag[8];
	int back_time_flag[8];
};


class DssdNormalizeEvent : public Event {
public:

	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	/// @param[in] prefix prefix of variables for friend tree
	///
	virtual void SetupInput(
		TTree *tree,
		const std::string &prefix = ""
	) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	unsigned short front_hit;
	unsigned short back_hit;
	unsigned short front_strip[8];
	unsigned short back_strip[8];
	double front_energy[8];
	double back_energy[8];
};


class DssdMergeEvent : public Event {
public:
	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	/// @param[in] prefix prefix of variables for friend tree
	///
	virtual void SetupInput(
		TTree *tree,
		const std::string &prefix = ""
	) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	unsigned short hit;
	unsigned int case_tag;
	double x[4];
	double y[4];
	double z[4];
	double energy[4];
	int time_flag[4];
};


class AdssdMergeEvent : public Event {
public:
	/// @brief setup branches of input tree
	/// @param[in] tree pointer to input tree
	/// @param[in] prefix prefix of variables for friend tree
	///
	virtual void SetupInput(
		TTree *tree,
		const std::string &prefix = ""
	) override;


	/// @brief setup branches of output tree
	/// @param[in] tree pointer to output tree
	///
	virtual void SetupOutput(TTree *tree) override;


	unsigned short hit;
	double radius[4];
	double theta[4];
	double phi[4];
	double energy[4];
	long long decode_entry[4];
};

}

#endif		// __DSSD_EVENT_H__