#ifndef __CHANNEL_H__
#define __CHANNEL_H__

#include <string>
#include <vector>

namespace ribll {

class Channel {
public:

	Channel(
		unsigned int run,
		const std::vector<std::string> &particles,
		const std::string &recoil = ""
	);

	virtual int Coincide();
protected:
	unsigned int run_;
	std::vector<std::string> fragments_;
	std::vector<unsigned short> charges_;
	std::vector<unsigned short> masses_;
	std::string recoil_;
	unsigned short recoil_mass_;
	unsigned short recoil_charge_;
};


class T0TAFChannel : public Channel {
public:
	T0TAFChannel(
		unsigned int run,
		const std::vector<std::string> &particles,
		const std::string &recoil = ""
	);

	virtual int Coincide() override;

};

class T0Channel : public Channel {
public:
	T0Channel(
		unsigned int run,
		const std::vector<std::string> &particles,
		const std::string &recoil = ""
	);

	virtual int Coincide() override;
};

}	// namespace ribll

#endif // __CHANNEL_H__