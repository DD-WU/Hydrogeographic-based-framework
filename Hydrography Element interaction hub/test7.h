#pragma once
#include"Feature_Interaction_Hub.h"
#include "Config.h"
#include"Boundary.h"
#include"River2D.h"
#include"land.h"
#include"Dam.h"
class test7:public Abstract_Feature_Interaction_Hub
{
public:
	test7(const char* config_name, const char* sheet_name);

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
	Boundary* boundary;
	River2D* river;
	Region* land;
	Dam* dam;
};

