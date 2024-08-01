# 模拟
## 编译好的程序在哪？
如果是按照前面的步骤用 cmake 生成的，就在 build/standalone/simulate/ 目录下。
## 需要运行哪些程序？
`generate.cpp` 生成数据，`t0_detect.cpp` 模拟了 T0 探测器阵列探测粒子；`taf_detect.cpp` 模拟了 TAF 探测器阵列探测粒子；`vtof_detect.cpp` 模拟了塑闪探测束流粒子，用于和 PPAC 配合使用（也可以不用）；`ppac_detect.cpp` 模拟了 PPAC 探测束流粒子。运行上述程序后，可以使用自己的程序分析。
也可以使用 `detect.cpp` 和 `simulate_rebuild.cpp` 进行快速分析，直达 Q 值谱，但是并不能用于检验自己的分析程序。
## 数据放在哪里？
自行修改 `include/defs.h` 中的 `kGenerateDataPath` 变量，该变量指示了数据的路径。在该路径下需要新建一些文件夹，以保存不同层级的数据。模拟至少应该包括 simulate、fundamental、normalize、merge、show、energy_calculate 文件夹。其中，energy_calculate 中的数据也要全盘复制，或者自行修改 build/standalone/calculate 下的程序并重新生成能损数据，还要提前导入 Lise++ 的数据。
