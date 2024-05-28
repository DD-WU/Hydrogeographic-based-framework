#pragma once
#include"Config.h"
#include <sstream>
#include"BasicFeature.h"
struct Node_Link
{
	int node;
	int link;
};
class Pipe_network : public BasicFeature
{
public:
	//flow functions
	Water_API Pipe_network(Config_pipe* parameter);
	Water_API void init()override;
	Water_API void update()override;
	Water_API void finalize()override;
	//inner functions
	Water_API void Load_pipe_upstreamlink_index();
	Water_API void Load_pipe_land_index();

	//data exchange functions
	Water_API Config_pipe* get_config();
	Water_API char* get_pipe_name(int i);
	Water_API int get_pipe_number();
	Water_API double get_pipe_elevation(int j);
	Water_API map<string, Node_Link> get_pipe_upstreamlink_index();
	Water_API map<string, int> get_pipe_land_index();
	Water_API double* get_overflow();
	Water_API void set_overflow(int i, double temp);
	Water_API double* get_inflow();
	Water_API void set_inflow(int i, double temp);
private:
	Config_pipe* parameter;
	char* pipe_name;
	double* p_val;//溢流数组
	double* q_in;//入流数组	
	map<string, Node_Link> pipe_upstreamlink_index;
	map<string, int> pipe_land_index;
};



