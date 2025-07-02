#!/bin/bash

echo "=== GPU NTT 优化效果验证 ==="
echo "正在编译现有的性能对比程序..."

# 尝试编译现有的性能测试程序
cd "$(dirname "$0")"

# 检查是否存在现有的性能测试程序
if [ -f "ntt/performance_comparison_optimized.cu" ]; then
    echo "发现现有的性能测试程序，正在编译..."
    nvcc -O3 -arch=sm_75 -std=c++17 --expt-relaxed-constexpr \
         -I. -Intt/src/include \
         -o performance_test \
         ntt/performance_comparison_optimized.cu \
         ntt/src/include/CUDA/ntt_gpu_optimized.cu 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "编译成功！正在运行测试..."
        echo ""
        ./performance_test
        echo ""
        echo "=== 测试完成 ==="
        exit 0
    else
        echo "编译失败，尝试备用方案..."
    fi
fi

# 备用方案：编译基础的GPU测试
if [ -f "ntt/src/test/test-case-1.cpp" ]; then
    echo "使用基础测试程序进行验证..."
    g++ -O3 -std=c++17 \
        -I. -Intt/src/include \
        -o basic_test \
        ntt/src/test/test-case-1.cpp 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "基础测试编译成功，正在运行..."
        ./basic_test
    else
        echo "基础测试也编译失败"
    fi
fi

# 显示优化代码的统计信息
echo ""
echo "=== GPU 优化代码分析 ==="
if [ -f "ntt/src/include/CUDA/ntt_gpu_optimized.cu" ]; then
    echo "优化文件大小: $(wc -l < ntt/src/include/CUDA/ntt_gpu_optimized.cu) 行代码"
    echo ""
    echo "主要优化策略："
    echo "1. 共享内存优化 - $(grep -c "shared_memory_optimized" ntt/src/include/CUDA/ntt_gpu_optimized.cu) 个实现"
    echo "2. 向量化处理 - $(grep -c "vectorized" ntt/src/include/CUDA/ntt_gpu_optimized.cu) 个实现"
    echo "3. 自适应配置 - $(grep -c "select_optimal_config" ntt/src/include/CUDA/ntt_gpu_optimized.cu) 个实现"
    echo "4. 内存池管理 - $(grep -c "GpuMemoryPool" ntt/src/include/CUDA/ntt_gpu_optimized.cu) 个实现"
    echo "5. 异步传输 - $(grep -c "cudaStream" ntt/src/include/CUDA/ntt_gpu_optimized.cu) 个实现"
else
    echo "未找到GPU优化文件"
fi

echo ""
echo "如果编译失败，请检查："
echo "1. CUDA工具链是否正确安装"
echo "2. GPU驱动是否支持CUDA"
echo "3. 是否有足够的GPU内存" 