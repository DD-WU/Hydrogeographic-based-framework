#pragma once
#include"Feature_Interaction_Hub.h"
#include"BasicFeature.h"
#include "Config.h"
#include"Boundary.h"
#include"River.h"
#include"land.h"
#include"Lake.h"
#include"Dam.h"
#include <iostream>

//只在这里试一下并行，其他方案没加(已测试，不好用，不加了)
class test2 :public Abstract_Feature_Interaction_Hub
{
public:
	test2(const char* config_name, const char* sheet_name);
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
	//std::vector<std::variant<std::shared_ptr<BasicFeature>>> objects;
	const char* config_name;
	const char* sheet_name;
	Config* config;
	Boundary* boundary;
	Region* land;
	River* river;
	Lake* lake;

};

