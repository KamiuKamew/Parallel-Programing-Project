# MPI NTT 性能测试脚本说明

本目录包含多个性能测试脚本，用于全面评估 MPI NTT 实现的性能特征。

## 脚本概览

### 1. `quick_test.sh` - 快速测试脚本

**用途**: 日常开发中快速验证功能和基本性能
**特点**:

- 运行时间短（约 1-2 分钟）
- 验证基本功能正确性
- 对比串行、2 进程、4 进程的基本性能

**使用方法**:

```bash
chmod +x scripts/quick_test.sh
./scripts/quick_test.sh
```

### 2. `scalability_test.sh` - 可扩展性测试脚本

**用途**: 深入分析强可扩展性和弱可扩展性
**特点**:

- 强可扩展性测试（固定问题规模，增加进程数）
- 弱可扩展性测试（按比例增加问题规模和进程数）
- 不同模数的可扩展性对比
- 混合并行（MPI + OpenMP）测试
- 生成 CSV 数据文件和分析报告

**使用方法**:

```bash
chmod +x scripts/scalability_test.sh
./scripts/scalability_test.sh
```

**输出文件**:

- `scalability_results/scalability_YYYYMMDD_HHMMSS.csv` - 详细数据
- `scalability_results/scalability_report_YYYYMMDD_HHMMSS.md` - 分析报告

### 3. `performance_test.sh` - 综合性能测试脚本

**用途**: 全面的性能测试套件
**特点**:

- 串行性能基准测试
- MPI 可扩展性测试
- 混合并行测试
- 通信开销分析
- 性能 profiling
- 自动生成详细报告

**使用方法**:

```bash
chmod +x scripts/performance_test.sh
./scripts/performance_test.sh
```

**输出文件**:

- `performance_results/performance_test_YYYYMMDD_HHMMSS.log` - 详细日志
- `performance_results/performance_report_YYYYMMDD_HHMMSS.md` - 性能报告

## 测试环境要求

### 必需软件

- **MPI 环境**: OpenMPI 或 MPICH
- **编译器**: GCC (支持 C++11)，mpic++
- **系统工具**: bc（用于浮点运算）

### 验证环境

```bash
# 检查MPI
mpirun --version
mpic++ --version

# 检查工具
which bc
```

### 目录结构要求

脚本需要在项目根目录运行，要求以下目录结构：

```
项目根目录/
├── ntt/                    # NTT源代码目录
│   ├── main.cc            # 主程序文件
│   └── src/               # 源代码子目录
├── scripts/               # 测试脚本目录
│   ├── quick_test.sh
│   ├── scalability_test.sh
│   ├── performance_test.sh
│   └── README.md
└── data/                  # 测试数据目录（可选）
    ├── input1.txt
    ├── input2.txt
    └── input3.txt
```

## 测试策略说明

### 1. 算法策略测试

- **块划分 vs 循环划分**: 验证我们选择块划分的合理性
- **流水线算法**: 分析流水线方法的适用性
- **任务分配优化**: 评估不同任务分配策略的效果

### 2. 编程方法测试

- **阻塞通信 vs 非阻塞通信**: 分析通信模式选择
- **集合通信 vs 点对点通信**: 评估通信效率
- **MPI 线程支持级别**: 测试不同线程支持级别的影响

### 3. 可扩展性分析

- **强可扩展性**: 固定问题规模，增加进程数
- **弱可扩展性**: 按比例增加问题规模和进程数
- **混合并行**: MPI 进程数 × OpenMP 线程数的组合优化

### 4. 性能特征分析

- **通信开销**: 分析不同进程数下的通信 vs 计算时间比
- **内存使用**: 评估不同配置下的内存占用
- **缓存性能**: 分析缓存友好性和局部性

## 数据解读指南

### 加速比 (Speedup)

```
加速比 = 串行时间 / 并行时间
```

- 理想值 = 进程数
- 实际值通常小于理想值（受通信开销影响）

### 效率 (Efficiency)

```
效率 = 加速比 / 进程数
```

- 理想值 = 1.0 (100%)
- 效率下降表明可扩展性瓶颈

### 关键指标

1. **最优进程数**: 效率开始显著下降的临界点
2. **通信开销**: 影响可扩展性的主要因素
3. **负载均衡**: 各进程工作量分配的均匀程度

## 故障排除

### 常见问题

1. **MPI 环境问题**

   ```bash
   # 症状：command not found: mpirun
   # 解决：安装MPI环境
   sudo apt-get install openmpi-bin openmpi-dev  # Ubuntu
   ```

2. **编译失败**

   ```bash
   # 症状：编译器找不到头文件
   # 解决：检查MPI开发包
   sudo apt-get install libopenmpi-dev
   ```

3. **权限问题**

   ```bash
   # 症状：Permission denied
   # 解决：添加执行权限
   chmod +x scripts/*.sh
   ```

4. **bc 命令缺失**
   ```bash
   # 症状：bc: command not found
   # 解决：安装计算器
   sudo apt-get install bc
   ```

### 性能异常诊断

1. **加速比异常低**

   - 检查通信开销是否过大
   - 验证负载均衡是否合理
   - 分析问题规模是否太小

2. **内存使用异常**

   - 检查内存分配策略
   - 验证数据结构设计
   - 分析内存泄漏可能性

3. **结果不一致**
   - 验证算法正确性
   - 检查数据竞争问题
   - 确认同步机制有效性

## 结果分析建议

1. **多次运行**: 每个配置至少运行 3-5 次，取平均值
2. **排除异常值**: 去除明显的异常数据点
3. **环境隔离**: 确保测试期间系统负载稳定
4. **结果验证**: 对比不同配置的计算结果正确性

## 扩展指南

### 添加新测试用例

1. 在对应脚本中添加新的测试函数
2. 更新 CSV 头部（如果需要新列）
3. 在 main 函数中调用新测试
4. 更新报告生成逻辑

### 自定义测试参数

修改脚本顶部的配置变量：

```bash
# 例如：修改测试的进程数范围
local process_counts=(1 2 4 8 16 32)

# 例如：修改测试的问题规模
local n_values=(65536 131072 262144)
```

这些脚本为 MPI NTT 性能分析提供了完整的测试框架，支持从快速验证到深度分析的各种需求。
