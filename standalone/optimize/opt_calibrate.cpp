#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TChain.h>

#include <ceres/ceres.h>
#include <glog/logging.h>

#include "include/event/t0_event.h"
#include "include/event/particle_event.h"
#include "include/calculator/range_energy_calculator.h"
#include "include/calculator/delta_energy_calculator.h"

using namespace ribll;

class CostFunctor {
public:
    CostFunctor(
        int layer,
        elc::DeltaEnergyCalculator *calculator,
        double channel1,
        double channel2
    )
    : layer_(layer)
    , calculator_(calculator)
    , channel1_(channel1)
    , channel2_(channel2) {}


    bool operator()(
        const double * const param1,
        const double * const param2,
        double *residual
    ) const {
        // double total_energy = 0.0;
        // double calibrate_energy[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        // for (int i = 0; i < layer_; ++i) {
        //     calibrate_energy[i] = param[2*i] + param[2*i+1]*channel_[i];
        //     total_energy += calibrate_energy[i];
        // }

        // double calculate_energy[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        // double range = calculator_->Range(total_energy);
        // for (int i = 0; i < layer_; ++i) {
        //     if (range < t0_thickness[i]) {
        //         calculate_energy[i] = calculator_->Energy(range);
        //         break;
        //     } else {
        //         double residual_range = range - t0_thickness[i];
        //         double residual_energy = calculator_->Energy(residual_range);
        //         calculate_energy[i] = total_energy - residual_energy;
        //         total_energy = residual_energy;
        //     }
        // }

        // for (int i = 0; i < 6; ++i) {
        //     residual[i] = calibrate_energy[i] - calculate_energy[i]; 
        // }

        // first layer calibrated energy
        double e1 = param1[0] + param1[1]*channel1_;
        // second layer calibrated energy
        double e2 = param2[0] + param2[1]*channel2_; 
        if (e1 < 0.0 || e2 < 0.0) {
            residual[0] = 1000;
            return true;
        }
        // 
        double e = calculator_->Energy(layer_, e1);
        // // total energys
        // double total_energy = e1 + e2;
        // // full range
        // double range = calculator_->Range(total_energy);
        // // second layer range
        // double residual_range = range - t0_thickness[layer_];
        // if (residual_range < 0.0) return false;
        // // second layer calculated energy
        // double e = calculator_->Energy(residual_range);
        // // first layer calculated energy
        // double de = total_energy - e;

        residual[0] = e - e2;

        return true;
    }

private:
    int layer_;
    elc::DeltaEnergyCalculator *calculator_;
    double channel1_;
    double channel2_;
};

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " run end_run\n";
        return -1;
    }
    int run = atoi(argv[1]);
    int end_run = atoi(argv[2]);

	google::InitGoogleLogging(argv[0]);

    // input T0 telescope chain
	TChain t0_chain("t0", "chain of T0 events");
	for (int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		t0_chain.AddFile(TString::Format(
			"%s%st0-telescope-ta-%04d.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			i
		));
	}
	// // input PPAC chain
	// TChain ppac_chain("ppac", "chain of ppac event");
	// for (int i = run; i <= end_run; ++i) {
	// 	if (i == 628) continue;
	// 	ppac_chain.AddFile(TString::Format(
	// 		"%s%sxppac-particle-ta-%04d.root/tree",
	// 		kGenerateDataPath,
	// 		kParticleDir,
	// 		i
	// 	));
	// }
	// // add friend
	// t0_chain.AddFriend(&ppac_chain, "ppac");
	// input T0 telescope event
	T0Event t0_event;
	// // input PPAC event
	// ParticleEvent ppac_event;
	// setup input branches
	t0_event.SetupInput(&t0_chain);
	// ppac_event.SetupInput(&t0_chain, "ppac.");

    // elc::RangeEnergyCalculator be_calculator("10Be", "Si");
    // elc::RangeEnergyCalculator he_calculator("4He", "Si");
    elc::DeltaEnergyCalculator be_calculator("t0", "10Be");
    elc::DeltaEnergyCalculator he_calculator("t0", "4He");

    ceres::Problem problem;

    double parameter[6][2] = {
        {0.0, 0.006},
        {0.0, 0.006},
        {0.0, 0.006},
        {0.0, 0.006},
        {0.0, 0.006},
        {0.0, 0.006}
    };

    // total number of entries
	long long entries = t0_chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Reading events   0%%");
	fflush(stdout);
	// loop to read events
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		t0_chain.GetEntry(entry);
		for (int i = 0; i < t0_event.num; ++i) {
            elc::DeltaEnergyCalculator *calculator = nullptr;
            if (t0_event.mass[i] == 10 && t0_event.charge[i] == 4) {
                calculator = &be_calculator;
            } else if (t0_event.mass[i] == 4 && t0_event.charge[i] == 2) {
                calculator = &he_calculator;
            } else {
                continue;
            }
            // double channel[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            // channel[0] = t0_event.energy[i][0];
            // channel[1] = t0_event.energy[i][1];
            // if (t0_event.layer[i] >= 2) channel[2] = t0_event.energy[i][2];
            // for (int j = 0; j < t0_event.layer[i]-2; ++j) {
            //     channel[j] = t0_event.ssd_energy[j];
            // }
            double channel1 = t0_event.layer[i] <= 3
                ? t0_event.energy[i][t0_event.layer[i]-1]
                : t0_event.ssd_energy[t0_event.layer[i]-4];
            double channel2 = t0_event.layer[i] <= 2
                ? t0_event.energy[i][t0_event.layer[i]]
                : t0_event.ssd_energy[t0_event.layer[i]-3];
            ceres::CostFunction *cost_func;
            cost_func = new ceres::NumericDiffCostFunction<
                CostFunctor, ceres::CENTRAL, 1, 2, 2
            > (
                new CostFunctor(
                    t0_event.layer[i]-1,
                    calculator,
                    channel1,
                    channel2
                )
            );
            problem.AddResidualBlock(cost_func, nullptr, parameter[t0_event.layer[i]-1], parameter[t0_event.layer[i]]);
        }
	}
	// show finish
	printf("\b\b\b\b100%%\n");

    ceres::Solver::Options options;
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";

    for (int i = 0; i < 6; ++i) {
        std::cout << parameter[i][0] << ", " << parameter[i][1] << "\n";
    }

    return 0;
}