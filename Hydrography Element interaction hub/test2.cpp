#include "test2.h"
#include <omp.h>
//#include <omp.h>


test2::test2(const char* config_name, const char* sheet_name)
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

void test2::beforeInit()
{
	//omp_set_num_threads(1); // 设置线程数量为 1
	//omp_set_dynamic(0); // 关闭动态线程数调整
	config = new Config(config_name, sheet_name);
	boundary = new Boundary(config);
	land = new Region(config);
	river = new River(config);
	lake = new Lake(config);
	//objects.push_back(std::shared_ptr<BasicFeature>(land));
	//objects.push_back(std::shared_ptr<BasicFeature>(river));
}

void test2::init()
{
	config->init();
	boundary->init();
	land->init();
	river->init();
	lake->init();
}

void test2::afterInit()
{
}

void test2::beforeUpdate()
{
}
// 289990ms加了更慢了
//258485ms
// 373103ms
void test2::update()
{
	land->update();//update自带并行，比外层加并行嵌套更快
	river->update();
//   // 获取系统可用的线程数
//	int num_threads = omp_get_max_threads();
//
//	// 并行更新对象
//#pragma omp parallel for
//	for (int i = 0; i < objects.size(); ++i) {
//		auto& obj = objects[i];
//#pragma omp parallel num_threads(threads_per_object) // 在并行区域内设置线程数
//		{
//			std::get<std::shared_ptr<BasicFeature>>(obj)->update();
//		}
//	}
//#pragma omp parallel for
//	for (int i = 0; i < objects.size(); ++i) {
//		auto& obj = objects[i];
//		std::get<std::shared_ptr<BasicFeature>>(obj)->update();
//	}
	config->update_time();
}

void test2::afterUpdate()
{
}

void test2::beforeFinalize()
{
}

void test2::finalize()
{
	land->finalize();
	river->finalize();
	config->finalize();
}

void test2::afterFinalize()
{
}

