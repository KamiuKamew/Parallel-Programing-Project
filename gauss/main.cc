#include <vector>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <sys/time.h>
#include <omp.h>
// 自行添加需要的头文件

int main(int argc, char *argv[])
{
    auto Start = std::chrono::high_resolution_clock::now();
    // TODO : 完成你的高斯消元
    // 如果需要文件输入输出, 请使用相对路径并放在 files 文件夹下
    auto End = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
    std::cout<<"average latency  : "<<elapsed.count()<<" (us) "<<std::endl;
    return 0;
}
