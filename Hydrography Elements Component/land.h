#pragma once

#ifdef Water_API//宏重定义解决办法
#undef Water_API
#endif

#if defined _WIN32
#define Water_API __declspec(dllexport)
/* Calling convention, stdcall in windows, cdecl in the rest of the world */
#define CALLCONV __stdcall
#else
#define Water_API
#define CALLCONV
#endif


#define MAXSTRINGLEN 1024
#define MAXDIMS 6
#include <stddef.h>
#include<sstream>
#include "Config.h"
class Region
{
public:
	//Work flow functions
	Water_API Region(const char* config, const char* sheet);
	Water_API Region(Config* lisf);
	Water_API void init();
	Water_API void update();
	Water_API void finalize();

	// inner fuctions
	Water_API void LoadStart(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose);
	Water_API void Load_land_pipe_index(Fnames* Fnameptr, States* Statesptr);
	Water_API void LoadManningsn(Fnames* Fnameptr, Pars* Parptr, Arrays* Arrptr, int* verbose);
	Water_API void FloodplainQ(States*, Pars*, Solver*, Arrays*);
	Water_API double CalcFPQx(int i, int j, States*, Pars*, Solver*, Arrays*, double* TSptr);
	Water_API double CalcFPQy(int i, int j, States*, Pars*, Solver*, Arrays*, double* TSptr);
 
	// Data exchange functions
	Water_API void set_inflow_lisflood(int j, double value);
	Water_API int get_gird_x(int j);
	Water_API int get_gird_y(int j);
	Water_API void set_waterdepth(int row, int col, double t);
	Water_API double get_elevation(int row, int col);
	Water_API double get_waterdepth(int row, int col);
	Water_API double* get_BCVar(int j);
	Water_API double get_dx();
	Water_API map<string, int> get_land_pipe_index();

private:
	Config* parameter;
	map< string, int > land_pipe_index;
};

//double Previous_t;      // previous time channel was calculated
//int steadyCount;
//int tstep_counter;   // start at -1 so that in first run through we calculate river



