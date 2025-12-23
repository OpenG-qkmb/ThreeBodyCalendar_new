# 三体人的万年历 ThreeBodyCalendar

这是计算机程序设计课程大作业。该项目拟实现模拟任意给定参数的三体运动，并以某一恒星附近的一颗行星为研究对象生成该行星的万年历。

This project is intended for submission as a project for the Computer Programming course. It aims to simulate the three-body motion of any given parameters, generating a perpetual calendar for a planet orbiting a specific star.


目录
- .gitattributes  
- .gitignore  
- CMakeLists.txt  
- CMakePresets.json  
- LICENSE.txt  
- README.md (此文件)  
- include/ (头文件目录)  
- src/ (实现代码)

项目目标
- 提供一个可配置的三体数值积分与后处理工具，用以：  
  - 判断给定初始条件下行星是否“稳定绕行”某恒星；  
  - 估计标准年（周期）并在长期演化中划分“恒纪元 / 乱纪元”；  
  - 导出粗略的“万年历”  

主要功能
- 数值积分模块（支持不同积分器的选择）  
- 即时与滑动窗口稳定性指标（能量和逃逸速度、Hill 判据等）  
- 标准年识别（尝试积分并采样，通过能量计算半长轴推知周期）  
- 纪元划分：通过瞬时判据周期性采样检测判断“恒纪元 / 乱纪元”  
- 输出：当前状态，主导恒星及纪元类型  
- 可视化前端：使用[EGE库](https://xege.org/)实现“示意图”

依赖
- CMake 3.15+  
- C++20 编译器  
- [Easy Graphics Engine](https://xege.org/)

算法与数值方法概述
- 引力模型：牛顿万有引力定律  
- 积分方法：
  - Euler法  
  - 速度Verlet法（默认）  
  - 四阶 Runge-Kutta 法  
- 稳定性判据：
  - 绑定判据：瞬时逃逸速度判断  
  - Hill 指标：计算行星与参考体的 Hill 半径并比较最近扰动天体距离  
  - ……
- 纪元划分逻辑：
  - 瞬时判据
  - 周期性采样

许可证
- 仓库内包含 LICENSE.txt，请参阅该文件获取具体许可证类型与使用条款。