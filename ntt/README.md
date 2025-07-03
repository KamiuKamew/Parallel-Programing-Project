# GPU NTT 并行化项目

## 项目简介

本项目实现了 Montgomery 规约 NTT 算法的 GPU 并行化优化，支持多尺度数据处理(n=4~512K)。

## 目录结构

```
ntt/
├── README.md                    # 项目说明文档
├── Makefile                     # 主编译文件
├── main.cc                      # 主程序入口
├── src/                         # 源代码目录
│   ├── cuda/                    # CUDA相关源码
│   │   ├── gpu_ntt_common.h     # GPU通用函数和数据结构
│   │   ├── gpu_ntt_v1.cu        # v1版本：预计算旋转因子优化
│   │   ├── gpu_ntt_v2.cu        # v2版本：共享内存优化
│   │   ├── gpu_ntt_v3.cu        # v3版本：多层共享内存优化
│   │   └── memory_wrapper.h     # RAII内存管理包装器
│   ├── include/                 # 头文件目录
│   │   └── ntt.h               # NTT算法接口声明
│   └── test/                    # 测试程序目录
│       ├── performance_test.cc  # 性能基准测试
│       └── correctness_test.cc  # 正确性验证测试
├── build/                       # 编译输出目录
├── docs/                        # 文档目录
│   ├── TODO                     # 任务规划文档
│   ├── optimization_summary.md  # 优化总结
│   └── performance_results.md   # 性能测试结果
└── scripts/                     # 构建和测试脚本
    ├── build.sh                 # 构建脚本
    └── test.sh                  # 测试脚本
```

## 编译和运行

```bash
# 编译主程序
make

# 运行性能测试
make test

# 清理编译产物
make clean
```

## 版本说明

- **v1**: 预计算旋转因子优化版本
- **v2**: 共享内存优化版本（推荐使用）
- **v3**: 多层共享内存实验版本

## 性能特点

- v2 版本在 n=131072 规模达到 78-86 微秒
- 相比串行版本提升约 12.5 倍
- 支持多尺度数据处理，算法正确性 100%验证
