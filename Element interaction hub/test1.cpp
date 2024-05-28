#include "test1.h"

test1::test1(const char* config_name, const char* sheet_name)
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

void test1::beforeInit()
{
	config = new Config(config_name, sheet_name);
	boundary = new Boundary(config);
	land = new Region(config);
	river = new River(config);
	lake = new Lake(config);
	dam = new Dam(config);
}

void test1::init()
{
	config->init();
	boundary->init();
	land->init();
	river->init();
	lake->init();
	dam->init();
}

void test1::afterInit()
{
}

void test1::beforeUpdate()
{
}

void test1::update()
{
	land->update();
	river->update();
	dam->update();
	config->update_time();
}

void test1::afterUpdate()
{
}

void test1::beforeFinalize()
{
}

void test1::finalize()
{
	land->finalize();
	river->finalize();
	config->finalize();
}

void test1::afterFinalize()
{
}
