#pragma once

#include <string>

#include "include/config.h"
#include "include/event/ingot/dssd_event.h"
#include "include/event/t0/dssd_match_event.h"

namespace brill {

struct DssdNormalizeParameters {
	int front_strips = 0;
	int back_strips = 0;
	double front_p0[kMaxStrips] = {0.0};
	double front_p1[kMaxStrips] = {1.0};
	double front_p2[kMaxStrips] = {0.0};
	double back_p0[kMaxStrips] = {0.0};
	double back_p1[kMaxStrips] = {1.0};
	double back_p2[kMaxStrips] = {0.0};
};

inline double NormEnergy(
	const brill::DssdNormalizeParameters &parameters,
	const int side,
	const int strip,
	const double raw_energy
) {
	if (side == 0) {
		return
			parameters.front_p0[strip]
			+ parameters.front_p1[strip] * raw_energy;
	}
	return
		parameters.back_p0[strip]
		+ parameters.back_p1[strip] * raw_energy;
}

int WriteDssdNormalizeParameters(
	const std::string &front_path,
	const std::string &back_path,
	const DssdNormalizeParameters &parameters
);

int ReadDssdNormalizeParameters(
	const std::string &front_path,
	const std::string &back_path,
	DssdNormalizeParameters &parameters
);

void ApplyDssdNormalize(
	const DssdEvent &input,
	const DssdNormalizeParameters &parameters,
	DssdEvent &output
);

void MatchDssdEvent(
	const DssdEvent &input,
	const SiliconDetectorConfig &detector,
	DssdMatchEvent &output
);

} // namespace brill
