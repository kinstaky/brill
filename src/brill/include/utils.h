#pragma once

#include <algorithm>
#include <string>

#include "include/config.h"

namespace brill {

inline std::string JoinPath(const std::string &left, const std::string &right) {
	if (left.empty()) return right;
	if (right.empty()) return left;
	if (left.back() == '/') return left + right;
	return left + "/" + right;
}

inline std::string TriggerInfix(const std::string &trigger) {
	return trigger == "main" ? "" : (trigger + "_");
}

inline bool IsJumpRun(const AppConfig &config, int run) {
	return std::find(config.jump_run.begin(), config.jump_run.end(), run) != config.jump_run.end();
}

} // namespace brill
