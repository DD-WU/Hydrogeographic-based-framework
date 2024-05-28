#include"test1.h"
#include"test2.h"
#include"test3.h"
#include"test4.h"
#include"test5.h"
#include"test6.h"
#include"test7.h"
#include"test8.h"
#include <chrono>
using namespace std::chrono;
int main()
{
	//////case A/case 1
	//new test1("F:\\SWMM_lisflood\\simulation\\parameter.xls", "test1");
	////case B
	//new test2("F:\\SWMM_lisflood\\simulation\\parameter.xls", "test2");
	////case C
	//new test3("F:\\SWMM_lisflood\\simulation\\parameter.xls", "test3");
	////case D
	//new test4("F:\\SWMM_lisflood\\simulation\\parameter.xls", "test4");
	////case E/case 4
	//new test5("F:\\SWMM_lisflood\\simulation\\parameter.xls", "test5");
	////case 3
	//new test6("F:\\SWMM_lisflood\\simulation\\parameter.xls", "test6");
	//case 2
	new test7("F:\\SWMM_lisflood\\simulation\\parameter.xls", "test7");
}

//
//auto start = std::chrono::high_resolution_clock::now(); // 获取开始时间点
//auto end = std::chrono::high_resolution_clock::now(); // 获取结束时间点
//// 计算代码运行时间并输出结果
//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//std::cout << "代码运行时间: " << duration.count() << " 毫秒" << std::endl;