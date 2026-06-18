#include "include/event/particle_event.h"

namespace brill {

void SetupInput(TTree *tree, ParticleEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "num").c_str(), &event.num);
	tree->SetBranchAddress((prefix + "charge").c_str(), event.charge);
	tree->SetBranchAddress((prefix + "mass").c_str(), event.mass);
	tree->SetBranchAddress((prefix + "energy").c_str(), event.energy);
	tree->SetBranchAddress((prefix + "time").c_str(), event.time);
	tree->SetBranchAddress((prefix + "x").c_str(), event.x);
	tree->SetBranchAddress((prefix + "y").c_str(), event.y);
	tree->SetBranchAddress((prefix + "z").c_str(), event.z);
	tree->SetBranchAddress((prefix + "px").c_str(), event.px);
	tree->SetBranchAddress((prefix + "py").c_str(), event.py);
	tree->SetBranchAddress((prefix + "pz").c_str(), event.pz);
	tree->SetBranchAddress((prefix + "stop").c_str(), event.stop);
	tree->SetBranchAddress((prefix + "last").c_str(), event.last);
}

void SetupOutput(TTree *tree, const ParticleEvent &event) {
	ParticleEvent &output = const_cast<ParticleEvent&>(event);
	tree->Branch("num", &output.num, "num/I");
	tree->Branch("charge", output.charge, "Z[num]/I");
	tree->Branch("mass", output.mass, "A[num]/I");
	tree->Branch("energy", output.energy, "e[num]/D");
	tree->Branch("time", output.time, "t[num]/D");
	tree->Branch("x", output.x, "x[num]/D");
	tree->Branch("y", output.y, "y[num]/D");
	tree->Branch("z", output.z, "z[num]/D");
	tree->Branch("px", output.px, "px[num]/D");
	tree->Branch("py", output.py, "py[num]/D");
	tree->Branch("pz", output.pz, "pz[num]/D");
	tree->Branch("stop", output.stop, "stop[num]/O");
	tree->Branch("last", output.last, "last[num]/I");
}

void Reset(ParticleEvent &event) {
	event.num = 0;
	for (int i = 0; i < 8; ++i) {
		event.charge[i] = 0;
		event.mass[i] = 0;
		event.energy[i] = 0.0;
		event.time[i] = 0.0;
		event.x[i] = 0.0;
		event.y[i] = 0.0;
		event.z[i] = 0.0;
		event.px[i] = 0.0;
		event.py[i] = 0.0;
		event.pz[i] = 0.0;
		event.stop[i] = false;
		event.last[i] = -1;
	}
}

} // namespace brill
