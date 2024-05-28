#include "Pipe_network.h"

Pipe_network::Pipe_network(Config_pipe* param)
{
    parameter = param;
    if (parameter)
    {
        cout << "   Creating pipe network module..." << endl;
        cout << "   " << endl;

    }
    else
    {
        cout << " config_pipe is not find" << endl;
        abort();
    }

}
void Pipe_network::init()
{
    swmm_start(true);
    pipe_name= parameter->get_pipe_name();
    int size = parameter->get_pipe_number();
    if (parameter->isRiverCouple()) {
        Load_pipe_upstreamlink_index();
    }
    if (parameter->isLandCouple()) {
        Load_pipe_land_index();
    }
    p_val = new double[size] {0};
    q_in = new double[size] {0};
}

void Pipe_network::update()
{
    double t= parameter->get_Elapsed_Time();
    swmm_step(&t);
}

void Pipe_network::finalize()
{
    swmm_end();
    swmm_close();
}

void Pipe_network::Load_pipe_upstreamlink_index()
{
    ifstream file(parameter->get_river_pipe_index_name());
    if (!file.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return;
    }
    std::string line;
    int num=parameter->get_pipe_number();
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;

        // 使用空格分隔键和值
        if (iss >> key >> value) {
            // 在这里你可以对key和value进行处理，比如存储到map中

            for (int i = 0; i < num; i++)
            {
                if (strcmp(value.c_str(), pipe_name + i * 80) == 0) {
                    Node_Link temp;
                    temp.node = get_Node_Id_swmm((char*)value.c_str());
                    temp.link = get_upstream_Link_swmm((char*)value.c_str());
                    pipe_upstreamlink_index[value] = temp;
                }
            }
        }
    }
    file.close();
}

void Pipe_network::Load_pipe_land_index()
{
    ifstream file(parameter->get_land_pipe_index_name());
    if (!file.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return;
    }
    int num = parameter->get_pipe_number();
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;

        // 使用空格分隔键和值
        if (iss >> key >> value) {
            // 在这里你可以对key和value进行处理，比如存储到map中
            for (int i = 0; i < num; i++)
            {
                if (strcmp((pipe_name + i * 80), value.c_str()) == 0)
                {
                    pipe_land_index[value] = i;
                }
            }
        }
    }
    file.close();
}

Config_pipe* Pipe_network::get_config()
{
    return parameter;
}

char* Pipe_network::get_pipe_name(int i)
{
    return pipe_name+ i * 80;
}

int Pipe_network::get_pipe_number()
{
    return parameter->get_pipe_number();
}

double Pipe_network::get_pipe_elevation(int j)
{
    return parameter->get_pipe_elevation(j);
}

map<string, Node_Link> Pipe_network::get_pipe_upstreamlink_index()
{
    return pipe_upstreamlink_index;
}

map<string, int> Pipe_network::get_pipe_land_index()
{
    return pipe_land_index;
}

double* Pipe_network::get_overflow()
{
    return p_val;
}

void Pipe_network::set_overflow(int i, double temp)
{
    p_val[i] = temp;
}


double* Pipe_network::get_inflow()
{
    return q_in;
}

void Pipe_network::set_inflow(int i, double temp)
{
    q_in[i] = temp;
}

//double Pipe_network::get_current_time(double* dt)
//{
//
//    if (dt != nullptr)
//    {
//        *dt = parameter->get_current_time();
//    }
//    return *dt;
//}




