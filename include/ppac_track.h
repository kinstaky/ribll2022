#ifndef __PPAC_TRACK_H__
#define __PPAC_TRACK_H__

#include <cmath>
#include <iostream>

#include "include/defs.h"

namespace ribll {

constexpr double ppac_correct[2][3] = {
	{0, 2.23, 3.4},
	{0, -0.84, -1.78}
};

constexpr double all_ppac_correct[2][4] = {
	{0.0, 0.95, 2.23, 3.4},
	{0.0, -0.16, -0.84, -1.78}
};


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


/// @brief linear fit points
/// @param[in] x x position of points
/// @param[in] y y position of points
/// @param[in] n number of points
/// @param[out] k slope of the line
/// @param[out] b intercept of the line
/// @returns chi2 of fitting
///
double SimpleFit(
	const double *x,
	const double *y,
	const int n,
	double &k,
	double &b
) {
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
	// chi2
	double chi2 = 0.0;
	for (int i = 0; i < n; ++i) {
		double t = y[i] - k*x[i] - b;
		chi2 += t * t;
	}
	return chi2;
}


/// @brief single PPAC tracking, optimize deutron approximate method
/// @param[in] c14_kinetic kinetic energy of 14C
/// @param[in] d_kinetic kinetic energy of 2H
/// @param[in] ppac_z distance of PPAC from target
/// @param[in] c14_position 10Be position on T0
/// @param[in] d_position 2H position on TAFD
/// @param[in] d_position_v 2H position on TAFD other direction
/// @param[in] c15_position beam position on PPAC
/// @returns reaction point
///
double ApproximateTrack(
	double c14_kinetic,
	double d_kinetic,
	double ppac_z,
	double c14_position,
	double d_position,
	double d_position_v,
	double c15_position
) {
	// 14C parameter
	double a_c14 = sqrt(14.0 * c14_kinetic) / 100.0;
	// 2H parameter
	double a_d =
		sqrt(2.0 * d_kinetic)
		/ sqrt(
			135.0 * 135.0
			+ d_position * d_position
			+ d_position_v * d_position_v
		);
	// 15C parameter
	double a_c15 = -sqrt(15.0 * 430.0) / ppac_z;
	// calculate
	double numerator = a_c14 * c14_position;
	numerator += a_d * d_position;
	numerator += a_c15 * c15_position;
	double denominator = a_c14 + a_d + a_c15;
	return numerator / denominator;
}


/// @brief Track PPAC points and get the line in 15C(p, d)14C reaction
/// @param[in] flag point valid flag
/// @param[in] x point x position (actually is z)
/// @param[in] y point y position (actually is x or y)
/// @param[in] be_kinetic kinematic energy of 10Be
/// @param[in] he_kinetic kinematic energy of 4He
/// @param[in] d_kinetic kinematic energy of 2H
/// @param[in] ppac_z distance of PPAC from target
/// @param[in] c14_position 10Be position on T0
/// @param[in] d_position 2H position on TAFD
/// @param[in] d_position_v 2H position on TAFD other direction
/// @param[out] reaction_point reaction point
/// @returns used PPAC number
///
int TrackPpac(
	unsigned short flag,
	const double *x,
	const double *y,
	double c14_kinetic,
	double d_kinetic,
	double c14_position,
	double d_position,
	double d_position_v,
	double &reaction_point
) {
	int ppac_num = 0;
	// beam slope and reaction point
	double k, b;
	b = -1e5;
	if (flag == 0x1) {
		b = ApproximateTrack(
			c14_kinetic, d_kinetic, x[0],
			c14_position, d_position, d_position_v, y[0]
		);
		ppac_num = 1;
	} else if (flag == 0x2) {
		b = ApproximateTrack(
			c14_kinetic, d_kinetic, x[1],
			c14_position, d_position, d_position_v, y[1]
		);
		ppac_num = 1;
	} else if (flag == 0x4) {
		b = ApproximateTrack(
			c14_kinetic, d_kinetic, x[2],
			c14_position, d_position, d_position_v, y[2]
		);
		ppac_num = 1;
	} else if (flag == 0x3) {
		k = (y[0] - y[1]) / (x[0] - x[1]);
		b = y[0] - k * x[0];
		ppac_num = 2;
	} else if (flag == 0x5) {
		k = (y[0] - y[2]) / (x[0] - x[2]);
		b = y[0] - k * x[0];
		ppac_num = 2;
	} else if (flag == 0x6) {
		k = (y[1] - y[2]) / (x[1] - x[2]);
		b = y[1] - k * x[1];
		ppac_num = 2;
	} else if (flag == 0x7) {
		double sumx = 0.0;
		double sumy = 0.0;
		double sumxy = 0.0;
		double sumx2 = 0.0;
		for (int i = 0; i < 3; ++i) {
			sumx += x[i];
			sumy += y[i];
			sumxy += x[i] * y[i];
			sumx2 += x[i] * x[i];
		}
		k = (sumxy - sumx*sumy/3.0) / (sumx2 - sumx*sumx/3.0);
		b = (sumy - k*sumx) / 3.0;
		ppac_num = 3;
	}
	reaction_point = b;
	return ppac_num;
}


/// @brief Track PPAC points and get the line
/// @param[in] flag point valid flag
/// @param[in] x point x position (actually is z)
/// @param[in] y point y position (actually is x or y)
/// @param[out] k slope of the line
/// @param[out] b intercept of the line
/// @returns 0 if success, -1 otherwise
///
int TrackMultiplePpac(
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
		SimpleFit(x, y, 3, k, b);
	} else {
		return -1;
	}
	return 0;
}


/// @brief single PPAC tracking, optimize deutron relative approximate method
/// @param[in] be_kinetic kinematic energy of 10Be
/// @param[in] he_kinetic kinematic energy of 4He
/// @param[in] d_kinetic kinematic energy of 2H
/// @param[in] ppac_z distance of PPAC from target
/// @param[in] be_position 10Be position on T0
/// @param[in] he_position 4He position on T0
/// @param[in] d_position 2H position on TAFD
/// @param[in] d_position_v 2H position on TAFD other direction
/// @param[in] c_position beam position on PPAC
/// @returns reaction point
///
double DeutronRelativeApproximateTrack(
	double be_kinetic,
	double he_kinetic,
	double d_kinetic,
	double ppac_z,
	double be_position,
	double he_position,
	double d_position,
	double d_position_v,
	double c_position
) {
	// 10Be parameter
	double a_be = sqrt(
		(2.0 * mass_10be + be_kinetic) * be_kinetic
	) / 100.0;
	// 4He parameter
	double a_he = sqrt(
		(2.0 * mass_4he + he_kinetic) * he_kinetic
	) / 100.0;
	// 2H parameter
	double a_d = sqrt(
		(2.0 * mass_2h + d_kinetic) * d_kinetic
	) / sqrt(
		135.0 * 135.0
		+ d_position * d_position
		+ d_position_v * d_position_v
	);
	// 14C parameter
	double a_c = -sqrt(
		(2.0 * mass_14c + 385.0) * 385.0
	) / ppac_z;
	// calculate
	double numerator = a_be * be_position;
	numerator += a_he * he_position;
	numerator += a_d * d_position;
	numerator += a_c * c_position;
	double denominator = a_be + a_he + a_d + a_c;
	return numerator / denominator;
}


/// @brief Track PPAC points and get the line in 14C(d, d')10Be+4He reaction
/// @param[in] flag point valid flag
/// @param[in] x point x position (actually is z)
/// @param[in] y point y position (actually is x or y)
/// @param[in] be_kinetic kinematic energy of 10Be
/// @param[in] he_kinetic kinematic energy of 4He
/// @param[in] d_kinetic kinematic energy of 2H
/// @param[in] ppac_z distance of PPAC from target
/// @param[in] be_position 10Be position on T0
/// @param[in] he_position 4He position on T0
/// @param[in] d_position 2H position on TAFD
/// @param[in] d_position_v 2H position on TAFD other direction
/// @param[out] reaction_point reaction point
/// @returns used PPAC number
///
int TrackPpac(
	unsigned short flag,
	const double *x,
	const double *y,
	double be_kinetic,
	double he_kinetic,
	double d_kinetic,
	double be_position,
	double he_position,
	double d_position,
	double d_position_v,
	double &reaction_point
) {
	int ppac_num = 0;
	// beam slope and reaction point
	double k, b;
	b = -1e5;
	if (flag == 0x1) {
		b = DeutronRelativeApproximateTrack(
			be_kinetic, he_kinetic, d_kinetic, x[0],
			be_position, he_position, d_position, d_position_v, y[0]
		);
		ppac_num = 1;
	} else if (flag == 0x2) {
		b = DeutronRelativeApproximateTrack(
			be_kinetic, he_kinetic, d_kinetic, x[1],
			be_position, he_position, d_position, d_position_v, y[1]
		);
		ppac_num = 1;
	} else if (flag == 0x4) {
		b = DeutronRelativeApproximateTrack(
			be_kinetic, he_kinetic, d_kinetic, x[2],
			be_position, he_position, d_position, d_position_v, y[2]
		);
		ppac_num = 1;
	} else if (flag == 0x3) {
		k = (y[0] - y[1]) / (x[0] - x[1]);
		b = y[0] - k * x[0];
		ppac_num = 2;
	} else if (flag == 0x5) {
		k = (y[0] - y[2]) / (x[0] - x[2]);
		b = y[0] - k * x[0];
		ppac_num = 2;
	} else if (flag == 0x6) {
		k = (y[1] - y[2]) / (x[1] - x[2]);
		b = y[1] - k * x[1];
		ppac_num = 2;
	} else if (flag == 0x7) {
		SimpleFit(x, y, 3, k, b);
		ppac_num = 3;
	}
	reaction_point = b;
	return ppac_num;
}


} // ribll

#endif	// __PPAC_TRACK_H__