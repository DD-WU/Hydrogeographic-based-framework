#pragma once
#include"Feature_Interaction_Hub.h"
#include "Config.h"
#include"tool.h"
#include"Boundary.h"
#include"River.h"
#include"land.h"
#include"Lake.h"
#include"Dam.h"
#include"road.h"
#include"Pipe_network.h"
class test8:public Abstract_Feature_Interaction_Hub
{
public:
	test8(const char* config_name, const char* sheet_name);

	virtual void beforeInit()override;
	virtual void init()override;
	virtual void afterInit()override;
	virtual void beforeUpdate()override;
	virtual void update()override;
	virtual void afterUpdate()override;
	virtual void beforeFinalize()override;
	virtual void finalize()override;
	virtual void afterFinalize()override;

private:
	const char* config_name;
	const char* sheet_name;

	Config* config;
	Config_pipe* config_pipe;

	Boundary* boundary;

	Region* land;
	River* river;
	Lake* lake;
	Dam* dam;
	Road* road;
	Pipe_network* pipe;

};


