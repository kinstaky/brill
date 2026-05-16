#include "include/t0/dssd.h"

#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/forge/dssd_event.h"
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
	return trigger.empty() ? "" : (trigger + "_");
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("match_dssd", "Match normalized DSSD events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("w,window", "Override matching tolerance.", cxxopts::value<double>(), "window")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		)
		(
			"detector",
			"Detectors to match.",
			cxxopts::value<std::vector<std::string>>(),
			"detector"
		);
	options.parse_positional({"detector"});

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
	if (!result.count("detector")) {
		std::cerr << "Error: Missing detector positional arguments.\n";
		PrintUsage(options);
		return 1;
	}

	const std::set<std::string> allowed_detectors = {"t0d1", "t0d2", "t0d3", "t0d4"};
	std::vector<std::string> detectors = result["detector"].as<std::vector<std::string>>();
	for (const std::string &detector : detectors) {
		if (allowed_detectors.count(detector) == 0) {
			std::cerr << "Error: Unsupported detector " << detector << ".\n";
			return 1;
		}
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) {
		return 1;
	}
	if (result.count("trigger")) {
		config.trigger = result["trigger"].as<std::string>();
	}
	const int run = result["run"].as<int>();

	for (const std::string &detector_name : detectors) {
		const brill::SquareDetectorConfig *detector =
			brill::FindDetectorConfig(config, detector_name);
		if (!detector) {
			std::cerr << "Error: Detector " << detector_name << " is not found in config.\n";
			return 1;
		}
		double match_tolerance = detector->match_tolerance;
		if (result.count("window")) {
			match_tolerance = result["window"].as<double>();
		}
		brill::SquareDetectorConfig working_detector = *detector;
		working_detector.match_tolerance = match_tolerance;

		brill::DssdNormalizeParameters parameters;
		parameters.front_strips = detector->front_strips;
		parameters.back_strips = detector->back_strips;
		std::string normalize_dir = JoinPath(config.workspace, config.paths.normalize);
		std::string front_path = JoinPath(normalize_dir, detector_name + "_front.txt");
		std::string back_path = JoinPath(normalize_dir, detector_name + "_back.txt");
		if (brill::ReadDssdNormalizeParameters(front_path, back_path, parameters)) {
			return 1;
		}

		TString input_path = TString::Format(
			"%s/%s_%s%04d.root",
			JoinPath(config.workspace, config.paths.forge).c_str(),
			detector_name.c_str(),
			TriggerInfix(config.trigger).c_str(),
			run
		);
		TFile ipf(input_path, "read");
		TTree *ipt = (TTree*)ipf.Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from " << input_path << " failed.\n";
			return 1;
		}

		std::filesystem::create_directories(JoinPath(config.workspace, config.paths.match));
		TString output_path = TString::Format(
			"%s/%s_%s%04d.root",
			JoinPath(config.workspace, config.paths.match).c_str(),
			detector_name.c_str(),
			TriggerInfix(config.trigger).c_str(),
			run
		);
		TFile opf(output_path, "recreate");
		TTree opt("tree", "matched dssd");
		brill::DssdMatchEvent match_event;
		brill::SetupOutput(&opt, match_event);

		brill::DssdEvent raw_event;
		brill::DssdEvent normalized_event;
		brill::SetupInput(ipt, raw_event);
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);
			brill::ApplyDssdNormalize(raw_event, parameters, normalized_event);
			brill::MatchDssdEvent(normalized_event, working_detector, match_event);
			opt.Fill();
		}

		opf.cd();
		opt.Write();
		opf.Close();
		ipf.Close();
	}

	return 0;
}
