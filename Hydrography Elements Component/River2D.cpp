#include "River2D.h"

extern "C" void fvcom_init(char* cstring[]);
extern "C" void update_fvcom();
extern "C" void finalize_fvcom();

extern "C" void Get_StartTime(float* startTime);
extern "C" void Get_EndTime(float* startTime);
extern "C" void Get_current_time_fvcom(int* startTime);
extern "C" void Set_Obc(float* h, int i);
extern "C" void get_depth_fvcom(float* h, int i);// h 水位，i不规则格点索引
extern "C" void get_water_depth_fvcom(float* h, int i);// h 水位，i不规则格点索引
extern "C" void get_elevation_fvcom(float* h, int i);// h 水位，i不规则格点索引
River2D::River2D(Config* config)
{
	this->parameter = config;
}

void River2D::init()
{
	fvcom_init(parameter->para);
	if (parameter->Statesptr->river2d_couple == 1)
	{
		Load_river_land_index(parameter->Fnameptr, parameter->Statesptr);
	}
	Load_river_2D_Flag(parameter->Fnameptr, parameter->Statesptr, parameter->Parptr, parameter->Arrptr, parameter->verbose);
}

void River2D::update()
{
	update_fvcom();
}

void River2D::finalize()
{
	finalize_fvcom();
}

void River2D::Load_river_land_index(Fnames* Fnameptr, States* Statesptr)
{
	ifstream file(Fnameptr->land_river_index_name);
	if (!file.is_open()) {
		std::cerr << "无法打开文件" << std::endl;
		return;
	}
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		int i, j, fvcom_node;

		// 使用空格分隔键和值
		if (iss >> i >> j>>fvcom_node) {
			Raster_index temp;
			temp.i = i; 
			temp.j = j;
			river_land_index[temp]= fvcom_node;
		}
	}
	file.close();
}
void River2D::Load_river_2D_Flag(Fnames* Fnameptr, States* Statesptr, Pars* Parptr, Arrays* Arrptr, int* verbose) {
	FILE* fp;
	int i, j, gr;
	char dum[800];
	double no_data_value = -9999;

	fp = fopen(Fnameptr->River2dFlagName, "r");
	if (fp == NULL)
	{
		return;
	}
	Arrptr->RiverFlag = new int[Parptr->xsz * Parptr->ysz]();
	for (i = 0; i < 5; i++) fscanf(fp, "%s %s", dum, dum);
	fscanf(fp, "%s %lf", dum, &no_data_value);
	for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
	{
		fscanf(fp, "%d", Arrptr->RiverFlag + i + j * Parptr->xsz);
		// if no_data set depth to zero
		if ((int)Arrptr->RiverFlag[i + j * Parptr->xsz] == no_data_value) Arrptr->RiverFlag[i + j * Parptr->xsz] = 0;
	}
	fclose(fp);

	if (*verbose == ON) printf("Done.\n\n");

	return;
};
void River2D::get_StartTime(float* startTime)
{
	Get_StartTime(startTime);
}

void River2D::get_EndTime(float* startTime)
{
	Get_EndTime(startTime);
}

void River2D::get_current_time(int* startTime)
{
	Get_current_time_fvcom(startTime);
}

void River2D::Set_Obc(float* h, int i)
{
}

void River2D::get_depth(float* h, int i)
{
	get_depth_fvcom(h, i);
}

void River2D::get_water_depth(float* h, int i)
{
	get_water_depth_fvcom(h, i);
}

void River2D::get_elevation(float* h, int i)
{
	get_elevation_fvcom(h, i);
}

map< Raster_index, int> River2D::get_river_land_index()
{
	return river_land_index;
}
