#ifndef __CHANNEL_H__
#define __CHANNEL_H__

#include <string>
#include <vector>

namespace ribll {

class Channel {
public:

	Channel(unsigned int run);

	virtual int Coincide();
protected:
	unsigned int run_;
};


class Be8ToTwoAlphaChannel : public Channel {
public:

	Be8ToTwoAlphaChannel(unsigned int run);

	virtual int Coincide();
};


class C12ToThreeAlphaChannel : public Channel {
public:

	C12ToThreeAlphaChannel(unsigned int run);

	virtual int Coincide();
};


class C14ToBe10He4TwoBodyChannel : public Channel {
public:
	C14ToBe10He4TwoBodyChannel(unsigned int run);

	virtual int Coincide();
};


class C14ToBe10He4ThreeBodyChannel : public Channel {
public:
	C14ToBe10He4ThreeBodyChannel(unsigned int run);

	virtual int Coincide();
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