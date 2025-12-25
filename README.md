# 三体人的万年历 ThreeBodyCalendar

这是计算机程序设计课程大作业。该项目拟实现模拟任意给定参数的三体运动，并以某一恒星附近的一颗行星为研究对象生成该行星的万年历。

This project is intended for submission as a project for the Computer Programming course. It aims to simulate the three-body motion of any given parameters, generating a perpetual calendar for a planet orbiting a specific star.


### 目录

```text
.gitattributes  
.gitignore  
CMakeLists.txt  
CMakePresets.json  
LICENSE.txt  
README.md (此文件)  
include/ (头文件目录)  
src/ (实现代码)
```

### 项目目标
- 提供一个可配置的三体数值积分与后处理工具，用以：  
  - 判断给定初始条件下行星是否“稳定绕行”某恒星；  
  - 估计标准年（周期）并在长期演化中划分“恒纪元 / 乱纪元”；  
  - 导出粗略的“万年历”  

### 主要功能
- 数值积分模块（支持不同积分器的选择）  
- 即时与滑动窗口稳定性指标（能量和逃逸速度、Hill 判据等）  
- 标准年识别（尝试积分并采样，通过能量计算半长轴推知周期）  
- 纪元划分：通过瞬时判据周期性采样检测判断“恒纪元 / 乱纪元”  
- 输出：当前状态，主导恒星及纪元类型  
- 可视化前端：使用[ EGE 库](https://xege.org/)实现“示意图”

### 复现项目的依赖
- CMake 3.15+  
- C++20 编译器  
- [Easy Graphics Engine](https://xege.org/)
- 建议：在 Visual Studio 2022 / 2026 版本中构建CMake项目

### 算法与数值方法概述
- 引力模型：牛顿万有引力定律  
- 积分方法：
  - Euler 法  
  - 速度 Verlet 法（默认）  
  - 四阶 Runge-Kutta 法  
- 稳定性判据：
  - 绑定判据：瞬时逃逸速度判断  
  - Hill 指标：计算行星与参考体的 Hill 半径并比较最近扰动天体距离  
  - ……
- 纪元划分逻辑：
  - 瞬时判据
  - 周期性采样

### 运行说明
运行分为初始化和设置两个步骤。
- 初始化  
	- 手动添加行星  
	```text
	init m
	```
	该模式下，恒星和行星均可手动添加
	```text
	add [type = star / planet] [id (must be unique)] [mass] [position_x] [position_y] [position_z] [velocity_x] [velocity_y] [velocity_z]
	```
	对于恒星，还可采取随机添加的方法
	```text
	rand [id] [mass]
	```
	对于行星，还可采取指定绕某恒星添加的方法（圆轨道）
	```text
	add [type = planet] [id] [mass] orbit [id_star] [radius]
	```
	注意：上述天体添加时的单位经过换算。  
	添加恒星时，质量以**太阳质量**为单位，长度以**天文单位**为单位，速度以**天文单位每小时**为单位；添加行星时，质量以**地球质量**为单位，其余与添加恒星时相同。  
	随机添加恒星时，恒星将被添加到离坐标原点5天文单位内的空间，速度在任一坐标轴上的分量不超过0.0016天文单位每小时。  
	指定绕某恒星添加行星时，不一定能使行星稳定绕行指定的恒星。例如：当指定的绕行半径过大时，行星可能被邻近其他天体扰动而脱离计划的轨道。  
	添加天体时各天体的id应当唯一；为生成某行星的万年历，必须添加**至少一颗行星**，否则不允许结束初始化环节。  
	所有天体添加完成后，输入结束标志以进入下一步：
	```text
	end init
	```
	- 从文件读取行星
	```text
	init s
	input [filename]
	```
	文件内包含的内容与手动添加的格式相同。结束标志亦应包括在文件末尾。  
	示例（手动添加）：
	```text
	init m
	add star sun0 1 0 0 0 0 0 0
	rand sun1 0.1221
	rand sun2 3.5129e-06
	add planet earth 1 orbit sun0 1
	end init
	```
	对上述例子，如需改为从文件（例如`input.txt`）添加，则可如下处理：  
	输入：
	```text
	init s
	input input.txt
	```
	此时确保`input.txt`中包含必要内容，例如：
	```text
	add star sun0 1 0 0 0 0 0 0
	rand sun1 0.1221
	rand sun2 3.5129e-06
	add planet earth 1 orbit sun0 1
	end init
	```
- 指定模拟设置  
在该步骤下，应指定模拟运动时的各项设置。
	- 积分方法  
	```text
	method [method = euler/verlet/rk4]
	```
	该项未经指定时，默认的积分方法是速度 Verlet 法（`verlet`）。如需更换（如换为四阶 Runge-Kutta 法，即`rk4`），则输入：
	```text
	method rk4
	```
	- 模拟时长  
	```text
	timelen [timelen_options]
	```
	模拟时长可选择固定时长（单位为**小时**），例如：
	```text
	timelen 876600
	```
	这选择模拟876600小时，即100个地球年的天体运动。  
	模拟时长亦可选择无限，例如：
	```text
	timelen unlimited
	```
	则此时程序将持续模拟和输出，直到进程被手动结束。  
	如不指定此项，则默认模拟时长为100个地球年。
	- 输出方法  
	输出万年历的方法默认为输出到屏幕（控制台），这与输入下面的命令效果一致：
	```text
	print2screen
	```
	一般而言，如果需要收集输出结果或考虑保持较高的模拟和输出速度时，建议选择输出到文件的选项：
	```text
	print2file [filename]
	```
	提示：在选择模拟时长为无限时不要选择输出到文件，否则输出文件的大小可能超过预期。此外，如选择输出到屏幕，在输出完成后程序将立即终止，建议在命令行下运行该程序以查看屏幕上输出的万年历信息。
	- 指定待分析天体  
	需要指定一颗行星待分析天体，例如在添加一颗名为`earth`的行星后，可如下指定其为待分析天体：
	```text
	analyze [id_planet]
	```
	如果不指定此项，则默认在天体列表中查找到的第一个行星为待分析天体。
	- 可视化选项  
	默认不开启可视化。如需开启（会导致万年历输出速率降低），则输入：
	```text
	display [width] [height]
	```
	其中`[width]`和`[height]`分别为可视化窗口的宽和高，单位为**像素**。  
	亦可不予指定窗口大小：
	```text
	display
	```
	则此时窗口大小默认为**1024×768**。
- 启动模拟  
完成上述设置后，输入`start`或`begin`以启动模拟。

### 许可证
- 仓库内包含 LICENSE.txt，请参阅该文件获取具体许可证类型与使用条款。