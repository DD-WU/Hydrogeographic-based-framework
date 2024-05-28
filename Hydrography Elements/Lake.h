#pragma once
#include"Config.h"
#include"BasicFeature.h"
class Lake : public BasicFeature
{
public:
	// flow functions
	Water_API Lake(Config* param);
	Water_API void init() override;
	Water_API void update() override;
	Water_API void finalize() override;

	// inner functions
	Water_API void LoadLakeStartWaterDepth(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose);

private:
	Config* parameter;
	int startflag;
};

