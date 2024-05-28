#include "test4.h"

test4::test4(const char* config_name, const char* sheet_name)
{
	this->config_name = config_name;
	this->sheet_name = sheet_name;
	beforeInit();
	init();
	afterInit();
	while (config->get_current_time() < config->get_sim_time()||config_pipe->get_current_time()!=0)
	{
		beforeUpdate();
		update();
		afterUpdate();
	}
	beforeFinalize();
	finalize();
	afterFinalize();
}

void test4::beforeInit()
{
	config = new Config(config_name, sheet_name);
	config_pipe = new Config_pipe((char*)config_name, (char*)sheet_name);
	boundary = new Boundary(config);
	land = new Region(config);
	river = new River(config);
	lake = new Lake(config);
	dam = new Dam(config);
	road = new Road(config);
	pipe= new Pipe_network(config_pipe);
}

void test4::init()
{
	config->init();
	boundary->init();
	land->init();
	river->init();
	lake->init();
	dam->init();
	road->init();
	pipe->init();
}

void test4::afterInit()
{
}
//试了一下，并行没有map迭代器快
void test4::beforeUpdate()
{
	double t1 = config->get_current_time();
	double t2 = config_pipe->get_current_time();
	if (config->get_current_time()==config_pipe->get_current_time())
		//理论上讲应该每个要素都有各自的时间，但由于我们是基于LISFLOOD、SWMM写的，模型在构建的时候处理好了一些要素之间的时间匹配关系，
		//因此我们就不多此一举给每个要素设计独立的时间了、而是做一下SWMM要素与LISFLOOD中的要素的时间匹配
	{
		map< string, riverBoudaryIndex*>& river_pipe_index = river->get_river_pipe_index();
		map<string, int >& land_pipe_index = land->get_land_pipe_index();
		map<string, Node_Link>& pipe_upstreamlink_index = pipe->get_pipe_upstreamlink_index();
		map<string, int >& pipe_land_index = pipe->get_pipe_land_index();
		int size = pipe->get_pipe_number();
		double* p_val = pipe->get_overflow();
		double* q_in = pipe->get_inflow();

		// 并行化循环
#pragma omp parallel for
		for (int i = 0; i < size; i++) { 
			p_val[i] = 0; 
		}
		config_pipe->set_inflow(q_in);//给swmm赋入流q_in，并获取溢流p_val
		config_pipe->get_overflow(p_val);

		map< string, riverBoudaryIndex*>::iterator river_it;
//#pragma omp parallel for
		for (river_it = river_pipe_index.begin(); river_it != river_pipe_index.end(); ++river_it)
		{
			map< string, Node_Link>::iterator temp = pipe_upstreamlink_index.find(river_it->first);
			if (temp != pipe_upstreamlink_index.end()) {
				Node_Link t = temp->second;
				double elevation = river-> get_river_elevation(*river_it->second);
				double H = river->get_river_H(*river_it->second);
				double pipe_elevation =pipe->get_pipe_elevation(t.node);
				config_pipe->setOutletDepth(t.node, ( H+ elevation));
				double q = config_pipe-> get_Link_Flow(t.link);
				river->set_river_inflow(q, *river_it->second);
			}
			else
			{
				cout << "请检查河流和管网的映射关系，河流映射表里有管网映射表里没有的键,我没想好怎么处理，为了不报错，我就不给你执行了" << endl;
				abort();
			}
		}
		map< string, int>::iterator land_it;
//#pragma omp parallel for
		for (land_it = land_pipe_index.begin(); land_it != land_pipe_index.end(); ++land_it) {
			map< string, int>::iterator temp = pipe_land_index.find(land_it->first);
			if (temp != pipe_land_index.end()) {
				land->set_inflow_lisflood(land_it->second, p_val[temp->second] / (land->get_dx()));
			}
		}
		double q_j = 0;//溢流量
		for (int i = 0; i < size; i++) q_in[i] = 0;
//#pragma omp parallel for
		for (land_it = land_pipe_index.begin(); land_it != land_pipe_index.end(); ++land_it) {
			map< string, int>::iterator temp = pipe_land_index.find(land_it->first);
			if (temp != pipe_land_index.end()) {
				//get position
				int ix = land->get_gird_x(land_it->second);
				int iy = land->get_gird_y(land_it->second);
				double waterdepth = land->get_waterdepth(ix, iy);
				q_j = Tool::waterDepthToInflow(waterdepth, 0.6, 0.75); //  leaving LISFLOOD  
				q_in[land_it->second] = q_in[land_it->second] + q_j / (land->get_dx()); // into SWMM 
				if (p_val[temp->second] > 0)
				{
					q_in[land_it->second] = 0;
				}
				else
				{
					land->set_inflow_lisflood(land_it->second, -q_j / (land->get_dx()));/**/
				}
			}
		}
	}
}

void test4::update()
{
	pipe->update();
	land->update();
	river->update();
	dam->update();
	config->update_time();
}

void test4::afterUpdate()
{
	//do something
}

void test4::beforeFinalize()
{
	//do something
}

void test4::finalize()
{
	pipe->finalize();
	land->finalize();
	river->finalize();
	config->finalize();
}

void test4::afterFinalize()
{
	//do something
}
