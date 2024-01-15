#ifndef __OPTMIZE_UTILITIES_H__
#define __OPTMIZE_UTILITIES_H__

#include <cmath>
#include <iostream>

constexpr double u = 931.494;
constexpr double mass_1h = u * 1.0072764520;
constexpr double mass_2h = u * 2.0135531980;
constexpr double mass_4he = u * 4.0015060943;
constexpr double mass_10be = u * 10.0113403769;
constexpr double mass_14c = u * 13.9999505089;
constexpr double mass_15c = u * 15.0073077289;

constexpr double ppac_xz[3] = {-695.2, -454.2, -275.2};
constexpr double ppac_yz[3] = {-689.2, -448.2, -269.2};


/// @brief get momentum from kinetic energy
/// @param[in] mass mass of particle
/// @param[in] kinetic kinetic energy of particle
/// @returns momentum of particle
///
inline double MomentumFromKinetic(double mass, double kinetic) {
	return sqrt((2.0 * mass + kinetic) * kinetic);
}


/// @brief linear fit points
/// @param[in] x x position of points
/// @param[in] y y position of points
/// @param[out] k slope of the line
/// @param[out] b intercept of the line
///
void SimpleFit(const double *x, const double *y, double &k, double &b) {
	int n = 3;
	double sumx = 0.0;
	double sumy = 0.0;
	double sumxy = 0.0;
	double sumx2 = 0.0;
	for (int i = 0; i < n; ++i) {
		sumx += x[i];
		sumy += y[i];
		sumxy += x[i] * y[i];
		sumx2 += x[i] * x[i];
	}
	k = (sumxy - sumx*sumy/double(n)) / (sumx2 - sumx*sumx/double(n));
	b = (sumy - k*sumx) / double(n);
}


/// @brief Track PPAC points and get the line
/// @param[in] flag point valid flag
/// @param[in] x point x position (actually is z)
/// @param[in] y point y position (actually is x or y)
/// @param[out] k slope of the line
/// @param[out] b intercept of the line
///
void TrackPpac(
	unsigned short flag,
	const double *x,
	const double *y,
	double &k,
	double &b
) {
	if (flag == 0x3) {
		k = (y[0] - y[1]) / (x[0] - x[1]);
		b = y[0] - k * x[0];
	} else if (flag == 0x5) {
		k = (y[0] - y[2]) / (x[0] - x[2]);
		b = y[0] - k * x[0];
	} else if (flag == 0x6) {
		k = (y[1] - y[2]) / (x[1] - x[2]);
		b = y[1] - k * x[1];
	} else if (flag == 0x7) {
		SimpleFit(x, y, k, b);
	} else {
		std::cerr << "Error: Unknown flag " << flag << "\n";
	}
}


/// @brief caculate the x offset with asummed l0
/// @param[in] l0 assumed offset of PPAC 0
/// @param[out] l1 calculated offset of PPAC 1
/// @param[out] l2 calculated offsest of PPAC 2
///
inline void PpacOffsetX(double l0, double &l1, double &l2) {
	l1 = 0.6969 * l0 + 2.23;
	l2 = 0.4718 * l0 + 3.40;
}


/// @brief caculate the y offset with asummed l0
/// @param[in] l0 assumed offset of PPAC 0
/// @param[out] l1 calculated offset of PPAC 1
/// @param[out] l2 calculated offsest of PPAC 2
///
inline void PpacOffsetY(double l0, double &l1, double &l2) {
	l1 = 0.6946 * l0 - 0.84;
	l2 = 0.4678 * l0 - 1.78;
}

#endif	// __OPTMIZE_UTILITIES_H__