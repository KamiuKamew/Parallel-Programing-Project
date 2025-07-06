#!/bin/bash

echo "=== CUDA NTT性能对比测试 ==="
echo "测试时间: $(date)"
echo ""

# 编译所有版本
echo "正在编译所有CUDA版本..."

# 基础版本
make -f Makefile.lab5_optimized main_gpu
if [ $? -eq 0 ]; then
    echo "✅ 基础版本编译成功"
else
    echo "❌ 基础版本编译失败"
    exit 1
fi

# Stage1优化版本
make -f Makefile.lab5_optimized main_gpu_stage1
if [ $? -eq 0 ]; then
    echo "✅ Stage1优化版本编译成功"
else
    echo "❌ Stage1优化版本编译失败"
fi

# Stage2优化版本
make -f Makefile.lab5_optimized main_gpu_stage2
if [ $? -eq 0 ]; then
    echo "✅ Stage2优化版本编译成功"
else
    echo "❌ Stage2优化版本编译失败"
fi

# Stage3优化版本
make -f Makefile.lab5_optimized main_gpu_stage3
if [ $? -eq 0 ]; then
    echo "✅ Stage3优化版本编译成功"
else
    echo "❌ Stage3优化版本编译失败"
fi

echo ""
echo "=== 性能测试开始 ==="
echo ""

# 测试函数
test_version() {
    local version_name=$1
    local executable=$2
    local mode=$3
    
    echo "测试版本: $version_name"
    echo "模乘模式: $mode"
    echo "----------------------------------------"
    
    # 运行测试并捕获输出
    output=$(./$executable 2>&1)
    
    # 提取性能数据
    n4_time=$(echo "$output" | grep "n = 4" | awk '{print $NF}' | sed 's/(us)//')
    n131072_p1_time=$(echo "$output" | grep "n = 131072.*p = 7340033" | awk '{print $NF}' | sed 's/(us)//')
    n131072_p2_time=$(echo "$output" | grep "n = 131072.*p = 104857601" | awk '{print $NF}' | sed 's/(us)//')
    n131072_p3_time=$(echo "$output" | grep "n = 131072.*p = 469762049" | awk '{print $NF}' | sed 's/(us)//')
    
    # 检查正确性
    correctness=$(echo "$output" | grep "多项式乘法结果正确" | wc -l)
    
    echo "n=4, p=7340033: ${n4_time:-N/A} μs"
    echo "n=131072, p=7340033: ${n131072_p1_time:-N/A} μs"
    echo "n=131072, p=104857601: ${n131072_p2_time:-N/A} μs"
    echo "n=131072, p=469762049: ${n131072_p3_time:-N/A} μs"
    echo "正确性测试通过: $correctness/4"
    echo ""
}

# 测试各个版本
echo "1. 基础版本 (Montgomery模乘)"
test_version "基础版本" "main_gpu" "Montgomery"

echo "2. Stage1优化版本 (内存访问优化)"
test_version "Stage1优化" "main_gpu_stage1" "Montgomery"

echo "3. Stage2优化版本 (共享内存优化)"
test_version "Stage2优化" "main_gpu_stage2" "Montgomery"

echo "4. Stage3优化版本 (综合优化)"
test_version "Stage3优化" "main_gpu_stage3" "Montgomery"

echo "5. Lab5最终优化版本 (智能共享内存)"
test_version "Lab5最终优化" "main_gpu" "Montgomery"

echo "=== 测试完成 ==="
echo "测试时间: $(date)" 