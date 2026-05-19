#include "include/event/ingot/silicon_event.h"

namespace brill {

void SetupInput(TTree *tree, SiliconEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"valid").c_str(), &event.valid);
	tree->SetBranchAddress((prefix+"time").c_str(), &event.time);
	tree->SetBranchAddress((prefix+"energy").c_str(), &event.energy);
}

void SetupOutput(TTree *tree, SiliconEvent &event) {
	tree->Branch("valid", &event.valid, "v/O");
	tree->Branch("time", &event.time, "t/D");
	tree->Branch("energy", &event.energy, "e/I");
}

}
