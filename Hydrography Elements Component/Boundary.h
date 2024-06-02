#pragma once
#include"Config.h"
class Boundary
{
public:
	Water_API Boundary(const char* config, const char* sheet);
	Water_API Boundary(Config* lisf);
	Water_API void init();
	Water_API void LoadBCs(Fnames*, States*, Pars*, BoundCs*, Arrays*, int*);
	Water_API void LoadBCVar(Fnames*, States*, Pars*, BoundCs*, 
		ChannelSegmentType*, Arrays*, vector<ChannelSegmentType>*, int*);
	Water_API int get_PSNum();
	Water_API int* get_BCPSNum();

private:
	Config* parameter;
};

