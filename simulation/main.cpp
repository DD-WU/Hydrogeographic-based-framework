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
//auto start = std::chrono::high_resolution_clock::now(); // ��ȡ��ʼʱ���
//auto end = std::chrono::high_resolution_clock::now(); // ��ȡ����ʱ���
//// �����������ʱ�䲢������
//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//std::cout << "��������ʱ��: " << duration.count() << " ����" << std::endl;