#include <ctime>
#include <iostream>

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>

#include "include/defs.h"

using namespace ribll;

constexpr double ppac_xz[3] = {-695.2, -454.2, -275.2};
constexpr double ppac_yz[3] = {-689.2, -448.2, -269.2};
constexpr double offset[3] = {0.0, -1.23, 0.76};

void NormalFit(const double *x, double *y, double *dy) {
	constexpr int n = 3;
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
	double k = (sumxy - sumx*sumy/double(n)) / (sumx2 - sumx*sumx/double(n));
	double b = (sumy - k*sumx) / double(n);
	for (int i = 0; i < n; ++i) {
		dy[i] = b + k * x[i] - y[i];
	}
	// double chi2 = 0.0;
	// for (int i = 0; i < n; ++i) {
	// 	double t = y[i] - k*x[i] - b;
	// 	chi2 += t * t;
	// }
	return;
}


void FixFit(const double *x, double *y, double x0, double y0, double *dy) {
	constexpr int n = 2;
	double a[n];
	double b[n];
	double sumab = 0.0;
	double suma2 = 0.0;
	for (size_t i = 0; i < n; ++i) {
		a[i] = x[i] - x0;
		b[i] = y[i] - y0;
		sumab += a[i] * b[i];
		suma2 += a[i] * a[i];
	}
	double k = sumab / suma2;
	for (size_t i = 0; i < n;  ++i) {
		dy[i] = -b[i] + k * a[i];
	}
	return;
}

void EquationSet(double *p0, double *p1, double *result) {
	constexpr int n = 3;
	double cofactor[n*n];
	cofactor[0] = p1[4] * p1[8] - p1[5] * p1[7];
	cofactor[1] = p1[5] * p1[6] - p1[3] * p1[8];
	cofactor[2] = p1[3] * p1[7] - p1[4] * p1[6];
	cofactor[3] = p1[2] * p1[7] - p1[1] * p1[8];
	cofactor[4] = p1[0] * p1[8] - p1[2] * p1[6];
	cofactor[5] = p1[1] * p1[6] - p1[0] * p1[7];
	cofactor[6] = p1[1] * p1[5] - p1[2] * p1[4];
	cofactor[7] = p1[2] * p1[3] - p1[0] * p1[5];
	cofactor[8] = p1[0] * p1[4] - p1[1] * p1[3];
	double t = p1[0] * cofactor[0] + p1[1] * cofactor[1] + p1[2] * cofactor[2];
	for (int i = 0; i < n; ++i) {
		result[i] = 0.0;
		for (int j = 0; j < n; ++j) {
			result[i] += p0[j] * cofactor[j*3+i];
		}
		result[i] /= t;
	}
	std::cout << result[0] << " " << result[1] << " " << result[2] << "\n";
	return;
}


void GetNormalizeOffset(const double *x, double *dx, double *calc_offset) {
	constexpr int n = 3;
	// sum x
	double sumx = 0.0;
	// sum x^2
	double sumx2 = 0.0;
	for (int i = 0; i < n; ++i) {
		sumx += x[i];
		sumx2 += x[i] * x[i];
	}
	// average x
	double ax = sumx / double(n);
	// denominator of paramters
	double t = sumx2 - ax * sumx;

	// p1 parameter in matrix
	double p1[n*n];
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			p1[i*3+j] = (x[i]-ax) * (x[j]-ax) / t;
			p1[i*3+j] += i == j ? -2.0/3.0 : 1.0/3.0;
		}
	}

	double denominator = p1[4] * p1[8] - p1[5] * p1[7];
	// result of equation set
	// double result[n];
	// EquationSet(dx, p1, result);
	calc_offset[0] = (dx[1] * p1[8] - p1[5] * dx[2]) / denominator;
	calc_offset[1] = (p1[4] * dx[2] - dx[1] * p1[7]) / denominator;

	std::cout << "dx0 " << dx[0] << ", dx1 " << dx[1] << ", dx2 " << dx[2] << "\n"
		<< "l0 " << 0.0 << ", l1 " << calc_offset[0] << ", l2 " << calc_offset[1] << "\n";
	return;
}


int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%sppac-offset.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of dx for normal fitting
	TH1F hist_normal_dx[3] {
		TH1F("hndx0", "#Deltax0", 1000, -10, 10),
		TH1F("hndx1", "#Deltax1", 1000, -10, 10),
		TH1F("hndx2", "#Deltax2", 1000, -10, 10)
	};
	// histogram of dx for fix fitting
	TH1F hist_fix_dx[2] {
		TH1F("hfdx1", "#Deltax1", 1000, -10, 10),
		TH1F("hfdx2", "#Deltax2", 1000, -10, 10)
	};

	// initialize random number generator
	TRandom3 generator(std::time(nullptr));

	// delta x get from normal fitting
	double normal_dx[3];
	// delta x get from fix fitting
	double fix_dx[2];

	int events = 100'000;
	for (int i = 0; i < events; ++i) {
		// position x at z=-800
		double x_at_800 = generator.Gaus(0.0, 5.0);
		// position x at z=0
		double x_at_0 = generator.Gaus(0.0, 5.0);
		// reality k
		double real_k = (x_at_800 - x_at_0) / -800.0;
		// reality b
		double real_b = x_at_0;

		// reality x0,x1,x2
		double real_x[3];
		for (size_t j = 0; j < 3; ++j) {
			real_x[j] = real_b + real_k * ppac_xz[j];
		}
		// measured position
		double measure_x[3];
		for (size_t j = 0; j < 3; ++j) {
			measure_x[j] = real_x[j] + offset[j];
		}

		NormalFit(ppac_xz, measure_x, normal_dx);
		FixFit(ppac_xz+1, measure_x+1, ppac_xz[0], measure_x[0], fix_dx);

		for (size_t j = 0; j < 3; ++j) hist_normal_dx[j].Fill(normal_dx[j]);
		for (size_t j = 0; j < 2; ++j) hist_fix_dx[j].Fill(fix_dx[j]);
	}

std::cout << "fix dx " << fix_dx[0] << " " << fix_dx[1] << "\n";

	// calculated offset from normalize fitting
	double norm_offset[2];
	GetNormalizeOffset(ppac_xz, normal_dx, norm_offset);
	// calculated offset from fix fitting
	// double fit_offset[2];


	// // fit normal offset
	// for (int i = 0; i < 3; ++i) {
	// 	TF1 fx(TString::Format("fnx%d", i), "gaus", -5, 5);
	// 	fx.SetParameter(0, 500);
	// 	fx.SetParameter(1, 0.0);
	// 	fx.SetParameter(2, 1.0);
	// 	hist_normal_dx[i].Fit(&fx, "QR+");
	// 	std::cout << "Normal dx" << i << " " << fx.GetParameter(1) << "\n";
	// }
	// // fit fix offset
	// for (int i = 0; i < 2; ++i) {
	// 	TF1 fx(TString::Format("ffx%d", i), "gaus", -5, 5);
	// 	fx.SetParameter(0, 500);
	// 	fx.SetParameter(1, 0.0);
	// 	fx.SetParameter(2, 1.0);
	// 	hist_fix_dx[i].Fit(&fx, "QR+");
	// 	std::cout << "Fix dx" << i+1 << " " << fx.GetParameter(1) << "\n";
	// }

	// save histogram
	for (size_t i = 0; i < 3; ++i) hist_normal_dx[i].Write();
	for (size_t i = 0; i < 2; ++i) hist_fix_dx[i].Write();
	// close file
	opf.Close();
	return 0;
}