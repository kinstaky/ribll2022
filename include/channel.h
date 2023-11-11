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


class C15pdChannel : public Channel {
public:
	C15pdChannel(unsigned int run);

	virtual int Coincide();
};


class C14ToBe10He4H1ThreeBodyChannel : public Channel {
public:
	C14ToBe10He4H1ThreeBodyChannel(unsigned int run);

	virtual int Coincide();
};


class C14ToHe4He4He6Channel : public Channel {
public:

	C14ToHe4He4He6Channel(unsigned int run);

	virtual int Coincide();
};


int MergeTaf(unsigned int run);


}	// namespace ribll

#endif // __CHANNEL_H__