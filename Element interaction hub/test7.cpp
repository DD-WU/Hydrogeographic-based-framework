#include "test7.h"
test7::test7(const char* config_name, const char* sheet_name)
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

void test7::beforeInit()
{
	config = new Config(config_name, sheet_name);
	boundary = new Boundary(config);
	land = new Region(config);
	dam = new Dam(config);
	river = new River2D(config);
}

void test7::init()
{
	config->init();
	boundary->init();
	land->init();
	dam->init();
	river->init();
	
}

void test7::afterInit()
{
}

void test7::beforeUpdate()
{
	int ft = -1;
	river->get_current_time(&ft);
	double lt = config->get_current_time();
	if (lt == ft)
	{
		map< Raster_index, int> land_river_index = river->get_river_land_index();
		int index;
		map< Raster_index, int>::iterator it;
		int count=0;
		double tt = land->get_waterdepth(224, 7);
		for (auto it = land_river_index.begin(); it != land_river_index.end(); ++it) {
			// 'it->first' is the key (Raster_index), 'it->second' is the value (int)
			float t, t1; float elevation;
			elevation = land->get_elevation(it->first.i, it->first.j);
			river->get_water_depth(&t, it->second);
			river->get_depth(&t1, it->second);
			land->set_waterdepth(it->first.i, it->first.j, t1);
		}
		cout << ft << endl;
	}
}

void test7::update()//这里举了个不同时间步长的例子
{
	int ft = -1;
	river->get_current_time(&ft);
	double lt = config->get_current_time();
	if (ft<lt)
	{
		river->update();
	}
	else if(ft > lt)
	{
		double tt = land->get_waterdepth(224, 7);
		land->update();
		dam->update();
		double save_t = config->get_save_time();
		double cur_t = config->get_current_time();
		if (cur_t==0|| (int)cur_t % (int)save_t == 0)
		{
			map< Raster_index, int> land_river_index = river->get_river_land_index();
			for (std::map<Raster_index, int>::iterator it = land_river_index.begin(); it != land_river_index.end(); ++it) {
				// 'it->first' is the key (Raster_index), 'it->second' is the value (int)
				float t, t1;
				float elevation;
				elevation= land->get_elevation(it->first.i, it->first.j);
				river->get_water_depth(&t, it->second);
				river->get_depth(&t1, it->second);
				land->set_waterdepth(it->first.i, it->first.j, t1);

			}
		}
		config->update_time();
	}
	else
	{

		land->update();
		dam->update();
		river->update();
		double save_t = config->get_save_time();
		double cur_t = config->get_current_time();
		if (cur_t == 0 || (int)cur_t%  (int)save_t == 0)
		{
			map< Raster_index, int> land_river_index = river->get_river_land_index();
			for (std::map<Raster_index, int>::iterator it = land_river_index.begin(); it != land_river_index.end(); ++it) {
				// 'it->first' is the key (Raster_index), 'it->second' is the value (int)
				float t, t1;
				river->get_water_depth(&t, it->second);
				river->get_depth(&t1, it->second);
				land->set_waterdepth(it->first.i, it->first.j, t1);
			}
		}
		config->update_time();
	}

}

void test7::afterUpdate()
{
}

void test7::beforeFinalize()
{
}

void test7::finalize()
{
	land->finalize();
	config->finalize();
}

void test7::afterFinalize()
{
}
