#pragma once
#include"Config.h"
#include<sstream>
#include"BasicFeature.h"
#pragma comment(lib, "..//Debug//fvcom.lib") //这个lib文件我没想好该放哪里，用户要是有比较好的位置放这个，可以改一下
//我把这个类的配置类设置为config了，但应该另写一个config_river2d的，主要是偷懒了
struct Raster_index
{
	int i, j;
	// Define the less-than operator for sorting in the map
	bool operator<(const Raster_index& other) const {
		if (i < other.i) return true;
		if (i > other.i) return false;
		return j < other.j;
	}
};
class River2D : public BasicFeature
{
public:
	//flow functions
	Water_API River2D(Config* config);
	Water_API void init()override;
	Water_API void update()override;
	Water_API void finalize()override;
	// inner functions
	Water_API void Load_river_land_index(Fnames* Fnameptr, States* Statesptr);
	Water_API void Load_river_2D_Flag(Fnames* Fnameptr, 
		States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose);
	//data exchange functions
	Water_API void get_StartTime(float* startTime);
	Water_API void get_EndTime(float* startTime);
	Water_API void get_current_time(int* startTime);
	Water_API void Set_Obc(float* h, int i);
	Water_API  void get_depth(float* h, int i);// h 水位，i不规则格点索引
	Water_API  void get_water_depth(float* h, int i);// h 水位，i不规则格点索引
	Water_API  void get_elevation(float* h, int i);// h 水位，i不规则格点索引
	Water_API map< Raster_index, int> get_river_land_index();
private:

	Config* parameter;
	map< Raster_index, int> river_land_index;
};

