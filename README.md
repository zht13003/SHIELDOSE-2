# SHIELDOSE2
SHIELDOSE-2模型的C/C++版本，[源代码](http://smdc.sinp.msu.ru/index.py?nav=model-shieldose2)基于Fortran 77编写~~根本看不懂~~，具体说明文档可以参照[这里](https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1116.pdf)。

用[f2c软件](http://www.netlib.org/f2c/)转为了C~~完全看不懂~~，编译为[dose.dll](/TotalDose/dose.dll)链接库。由于源代码采用的是读取输入文件（目录下的[example.inp](/TotalDose/ConsoleApplication1/example.inp)）文件然后返回剂量值并输出结果文件（目录下的[example.OUT](/TotalDose/ConsoleApplication1/example.OUT)），写了一段创建输入文件的代码并编译为[dose_input.dll](/TotalDose/dose_input.dll)链接库，输入为等效铝屏蔽厚度、质子能谱，如有其他需求可另行编写。

具体示例可见[TotalDose.cpp](/TotalDose/ConsoleApplication1/TotalDose.cpp)
