#!/bin/bash
# GPU NTT项目构建脚本

echo "=== GPU NTT项目构建脚本 ==="

# 检查CUDA环境
if ! command -v nvcc &> /dev/null; then
    echo "❌ 错误: nvcc未找到，请安装CUDA toolkit"
    exit 1
fi

echo "✅ CUDA环境检查通过"
nvcc --version

# 创建构建目录
mkdir -p build

# 根据参数选择构建版本
VERSION=${1:-v2}

case $VERSION in
    v1)
        echo "🔨 构建v1版本（预计算旋转因子优化）"
        make v1
        ;;
    v2)
        echo "🔨 构建v2版本（共享内存优化，推荐）"
        make
        ;;
    v3)
        echo "🔨 构建v3版本（多层共享内存实验）"
        make v3
        ;;
    test)
        echo "🔨 构建性能测试程序"
        make test
        ;;
    clean)
        echo "🧹 清理构建产物"
        make clean
        ;;
    all)
        echo "🔨 构建所有版本"
        make v1 && make && make v3 && make test
        ;;
    *)
        echo "使用方法: $0 [v1|v2|v3|test|clean|all]"
        echo "  v1   - 构建v1版本（预计算旋转因子优化）"
        echo "  v2   - 构建v2版本（共享内存优化，推荐，默认）"
        echo "  v3   - 构建v3版本（多层共享内存实验）"
        echo "  test - 构建性能测试程序"
        echo "  clean- 清理构建产物"
        echo "  all  - 构建所有版本"
        exit 1
        ;;
esac

if [ $? -eq 0 ]; then
    echo "✅ 构建完成！"
    echo "📁 可执行文件位于 build/ 目录中"
else
    echo "❌ 构建失败！"
    exit 1
fi 