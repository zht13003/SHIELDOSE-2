//|外部文件:
//|	名称		      格式			大小				描述
//| elbrbas2.dat      txt           312K              数据文件
//|                                              
//| protbas2.dat      txt           325K              数据文件
//|                                              
//| vcf2c.lib        静态库         317K              编译所需的静态库

#include <iostream>
#include <Windows.h>
using namespace std;
int main()
{
    typedef void (*input)(double depth, double* energy, double* flux, int length);
    HMODULE hDll = LoadLibraryA("dose_input.dll");
    input p1 = input(GetProcAddress(hDll, "?input@doseModel@dose@@SAXNPAN0H@Z"));

    double depth = 5e-2;
    double e[5] = { 1,2,3,4,5 };
    double f[5] = { 1e9,1e8,1e7,1e6,1e5 };
    p1(depth, e, f, sizeof(e) / sizeof(e[0]));//创建输入文件example.inp
                                              //输入为等效铝屏蔽厚度（g/cm^2)、质子能量（MeV）、对应的质子积分通量（/cm^2*s)、数组长度

    typedef double (*getDose)();
    HMODULE hDll2 = LoadLibraryA("dose.dll");
    getDose p2 = getDose(GetProcAddress(hDll2, "_dose@0"));

    cout << p2();//读取输入文件，创建输出文件example.OUT并返回硅材料剂量值
}


