# ribll 2022 data analysis

A project to analysis data from experiemnt in ribll 1 in 2022 by PKU.

## 模拟

只需要模拟
```bash
cmake -B build -S . -DBUILD_ALL=OFF -DBUILD_SIMULATION=ON
cmake --build build -- -j8
```
第一句是用 cmake 生成各种 Makefile，第二句用于生成可执行文件。第一句只有第一次使用时需要，第二句每次更改完程序都要运行以重新编译，第二句等效于 `make -j8`。详细参见 [模拟](standalone/simulate/README.md)
