#include "test3.h"

test3::test3(const char* config_name, const char* sheet_name)
{
	this->config_name = config_name;
	this->sheet_name = sheet_name;
	beforeInit();
	init();
	afterInit();
	while (config->get_current_time() < config->get_sim_time())
	{
		beforeUpdate();
		update();
		afterUpdate();
	}
	beforeFinalize();
	finalize();
	afterFinalize();
}

void test3::beforeInit()
{
	config = new Config(config_name, sheet_name);
	boundary = new Boundary(config);
	land = new Region(config);
	river = new River(config);
	lake = new Lake(config);
	dam = new Dam(config);
	road = new Road(config);
}

void test3::init()
{
	config->init();
	boundary->init();
	land->init();
	river->init();
	lake->init();
	dam->init();
	road->init();
}

void test3::afterInit()
{
}

void test3::beforeUpdate()
{
}

void test3::update()
{
	land->update();
	river->update();
	dam->update();
	config->update_time();
}

void test3::afterUpdate()
{
}

void test3::beforeFinalize()
{
}

void test3::finalize()
{
	land->finalize();
	river->finalize();
	config->finalize();
}

void test3::afterFinalize()
{
}
