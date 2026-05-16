#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/t0/dssd_match_event.h"

namespace {

std::string JoinPath(const std::string &left, const std::string &right) {
	if (left.empty()) return right;
	if (right.empty()) return left;
	if (left.back() == '/') return left + right;
	return left + "/" + right;
}

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

std::string TriggerInfix(const std::string &trigger) {
	return trigger == "main" ? "" : (trigger + "_");
}

double AxisPitchX(const brill::SquareDetectorConfig &detector) {
	if (detector.name == "t0d2") {
		return detector.size_x_mm / double(detector.back_strips);
	}
	return detector.size_x_mm / double(detector.front_strips);
}

double AxisPitchY(const brill::SquareDetectorConfig &detector) {
	if (detector.name == "t0d2") {
		return detector.size_y_mm / double(detector.front_strips);
	}
	return detector.size_y_mm / double(detector.back_strips);
}

double RandomizeCoordinate(double value, double pitch, TRandom3 &generator) {
	return value + generator.Rndm() * pitch - 0.5 * pitch;
}

void FillOffsets(
	TH1F &hist_x,
	TH1F &hist_y,
	const brill::DssdMatchEvent &left,
	double left_pitch_x,
	double left_pitch_y,
	const brill::DssdMatchEvent &right,
	double right_pitch_x,
	double right_pitch_y,
	TRandom3 &generator
) {
	for (int i = 0; i < left.num; ++i) {
		double left_x = RandomizeCoordinate(left.x[i], left_pitch_x, generator);
		double left_y = RandomizeCoordinate(left.y[i], left_pitch_y, generator);
		for (int j = 0; j < right.num; ++j) {
			double right_x = RandomizeCoordinate(right.x[j], right_pitch_x, generator);
			double right_y = RandomizeCoordinate(right.y[j], right_pitch_y, generator);
			hist_x.Fill(right_x - left_x);
			hist_y.Fill(right_y - left_y);
		}
	}
}

void FitAndPrint(TH1F &histogram) {
	TF1 fit((std::string("fit_") + histogram.GetName()).c_str(), "gaus", -5.0, 5.0);
	fit.SetParameter(1, 0.0);
	fit.SetParameter(2, 1.0);
	if (histogram.GetEntries() > 0) {
		histogram.Fit(&fit, "QR+");
		std::cout
			<< histogram.GetName()
			<< ": mean = " << fit.GetParameter(1)
			<< ", sigma = " << fit.GetParameter(2)
			<< "\n";
	} else {
		std::cout
			<< histogram.GetName()
			<< ": empty histogram\n";
	}
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_t0_center", "Estimate T0 detector center offsets.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		);

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		PrintUsage(options);
		return 0;
	}
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		PrintUsage(options);
		return 1;
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) {
		return 1;
	}
	if (result.count("trigger")) {
		config.trigger = result["trigger"].as<std::string>();
	}

	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return 1;
	}

	const brill::SquareDetectorConfig *detector1 = brill::FindDetectorConfig(config, "t0d1");
	const brill::SquareDetectorConfig *detector2 = brill::FindDetectorConfig(config, "t0d2");
	const brill::SquareDetectorConfig *detector3 = brill::FindDetectorConfig(config, "t0d3");
	const brill::SquareDetectorConfig *detector4 = brill::FindDetectorConfig(config, "t0d4");
	if (!detector1 || !detector2 || !detector3 || !detector4) {
		std::cerr << "Error: Missing T0 detector config.\n";
		return 1;
	}

	const std::string match_dir = JoinPath(config.workspace, config.paths.match);
	const std::string trigger_infix = TriggerInfix(config.trigger);

	TChain chain1("tree");
	TChain chain2("tree");
	TChain chain3("tree");
	TChain chain4("tree");
	for (int current_run = run; current_run <= end_run; ++current_run) {
		chain1.Add(TString::Format(
			"%s/t0d1_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain2.Add(TString::Format(
			"%s/t0d2_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain3.Add(TString::Format(
			"%s/t0d3_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain4.Add(TString::Format(
			"%s/t0d4_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
	}
	chain1.AddFriend(&chain2, "d2");
	chain1.AddFriend(&chain3, "d3");
	chain1.AddFriend(&chain4, "d4");

	brill::DssdMatchEvent event1;
	brill::DssdMatchEvent event2;
	brill::DssdMatchEvent event3;
	brill::DssdMatchEvent event4;
	brill::SetupInput(&chain1, event1);
	brill::SetupInput(&chain1, event2, "d2.");
	brill::SetupInput(&chain1, event3, "d3.");
	brill::SetupInput(&chain1, event4, "d4.");

	TString output_path = TString::Format(
		"%s/t0_center_%s%04d_%04d.root",
		JoinPath(config.workspace, config.paths.estimate).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	TH1F d2d1_x_offset("d2d1_x_offset", "D2 x - D1 x", 1000, -10.0, 10.0);
	TH1F d2d1_y_offset("d2d1_y_offset", "D2 y - D1 y", 1000, -10.0, 10.0);
	TH1F d3d1_x_offset("d3d1_x_offset", "D3 x - D1 x", 1000, -10.0, 10.0);
	TH1F d3d1_y_offset("d3d1_y_offset", "D3 y - D1 y", 1000, -10.0, 10.0);
	TH1F d4d1_x_offset("d4d1_x_offset", "D4 x - D1 x", 1000, -10.0, 10.0);
	TH1F d4d1_y_offset("d4d1_y_offset", "D4 y - D1 y", 1000, -10.0, 10.0);
	TH1F d3d2_x_offset("d3d2_x_offset", "D3 x - D2 x", 1000, -10.0, 10.0);
	TH1F d3d2_y_offset("d3d2_y_offset", "D3 y - D2 y", 1000, -10.0, 10.0);
	TH1F d4d2_x_offset("d4d2_x_offset", "D4 x - D2 x", 1000, -10.0, 10.0);
	TH1F d4d2_y_offset("d4d2_y_offset", "D4 y - D2 y", 1000, -10.0, 10.0);
	TH1F d4d3_x_offset("d4d3_x_offset", "D4 x - D3 x", 1000, -10.0, 10.0);
	TH1F d4d3_y_offset("d4d3_y_offset", "D4 y - D3 y", 1000, -10.0, 10.0);

	const double d1_pitch_x = AxisPitchX(*detector1);
	const double d1_pitch_y = AxisPitchY(*detector1);
	const double d2_pitch_x = AxisPitchX(*detector2);
	const double d2_pitch_y = AxisPitchY(*detector2);
	const double d3_pitch_x = AxisPitchX(*detector3);
	const double d3_pitch_y = AxisPitchY(*detector3);
	const double d4_pitch_x = AxisPitchX(*detector4);
	const double d4_pitch_y = AxisPitchY(*detector4);

	const long long total = chain1.GetEntries();
	TRandom3 generator(total);
	long long last_percentage = -1;
	std::printf("Estimating T0 center   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);

		FillOffsets(
			d2d1_x_offset,
			d2d1_y_offset,
			event1,
			d1_pitch_x,
			d1_pitch_y,
			event2,
			d2_pitch_x,
			d2_pitch_y,
			generator
		);
		FillOffsets(
			d3d1_x_offset,
			d3d1_y_offset,
			event1,
			d1_pitch_x,
			d1_pitch_y,
			event3,
			d3_pitch_x,
			d3_pitch_y,
			generator
		);
		FillOffsets(
			d4d1_x_offset,
			d4d1_y_offset,
			event1,
			d1_pitch_x,
			d1_pitch_y,
			event4,
			d4_pitch_x,
			d4_pitch_y,
			generator
		);
		FillOffsets(
			d3d2_x_offset,
			d3d2_y_offset,
			event2,
			d2_pitch_x,
			d2_pitch_y,
			event3,
			d3_pitch_x,
			d3_pitch_y,
			generator
		);
		FillOffsets(
			d4d2_x_offset,
			d4d2_y_offset,
			event2,
			d2_pitch_x,
			d2_pitch_y,
			event4,
			d4_pitch_x,
			d4_pitch_y,
			generator
		);
		FillOffsets(
			d4d3_x_offset,
			d4d3_y_offset,
			event3,
			d3_pitch_x,
			d3_pitch_y,
			event4,
			d4_pitch_x,
			d4_pitch_y,
			generator
		);
	}
	std::printf("\b\b\b\b100%%\n");

	FitAndPrint(d2d1_x_offset);
	FitAndPrint(d2d1_y_offset);
	FitAndPrint(d3d1_x_offset);
	FitAndPrint(d3d1_y_offset);
	FitAndPrint(d4d1_x_offset);
	FitAndPrint(d4d1_y_offset);
	FitAndPrint(d3d2_x_offset);
	FitAndPrint(d3d2_y_offset);
	FitAndPrint(d4d2_x_offset);
	FitAndPrint(d4d2_y_offset);
	FitAndPrint(d4d3_x_offset);
	FitAndPrint(d4d3_y_offset);

	opf.cd();
	d2d1_x_offset.Write();
	d2d1_y_offset.Write();
	d3d1_x_offset.Write();
	d3d1_y_offset.Write();
	d4d1_x_offset.Write();
	d4d1_y_offset.Write();
	d3d2_x_offset.Write();
	d3d2_y_offset.Write();
	d4d2_x_offset.Write();
	d4d2_y_offset.Write();
	d4d3_x_offset.Write();
	d4d3_y_offset.Write();
	opf.Close();

	return 0;
}
