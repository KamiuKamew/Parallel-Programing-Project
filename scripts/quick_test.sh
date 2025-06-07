#!/bin/bash

# MPI NTT 快速测试脚本
# 用于快速验证功能和基本性能

set -e

# 颜色定义
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log() {
    echo -e "${GREEN}[$(date '+%H:%M:%S')] $1${NC}"
}

warn() {
    echo -e "${YELLOW}[$(date '+%H:%M:%S')] $1${NC}"
}

error() {
    echo -e "${RED}[$(date '+%H:%M:%S')] $1${NC}"
}

# 检查环境
if [ ! -d "ntt" ]; then
    error "请在项目根目录运行此脚本"
    exit 1
fi

if ! command -v mpirun &> /dev/null; then
    error "MPI环境未找到"
    exit 1
fi

log "开始快速测试..."

# 编译
cd ntt
log "编译代码..."
mpic++ -O3 -fopenmp -DUSE_MPI -o main_mpi main.cc || {
    error "编译失败"
    exit 1
}

# 快速功能测试
log "功能测试 - 小规模数据..."
echo "4 7340033" | ./main_mpi && log "✓ 小规模测试通过" || error "✗ 小规模测试失败"

# 性能对比测试
log "性能对比测试 - n=131072, 模数=104857601..."

# 串行基准（使用MPI但只启动1个进程）
log "测试串行版本（1进程MPI）..."
serial_output=$(time (echo "131072 104857601" | mpirun -np 1 ./main_mpi > /dev/null) 2>&1)
serial_time=$(echo "$serial_output" | grep real | awk '{print $2}')
log "串行时间: $serial_time"

# MPI 2进程
log "测试MPI 2进程..."
mpi2_output=$(time (echo "131072 104857601" | mpirun -np 2 ./main_mpi > /dev/null) 2>&1)
mpi2_time=$(echo "$mpi2_output" | grep real | awk '{print $2}')
log "MPI 2进程时间: $mpi2_time"

# MPI 4进程
log "测试MPI 4进程..."
mpi4_output=$(time (echo "131072 104857601" | mpirun -np 4 ./main_mpi > /dev/null) 2>&1)
mpi4_time=$(echo "$mpi4_output" | grep real | awk '{print $2}')
log "MPI 4进程时间: $mpi4_time"

cd ..

log "快速测试完成！"
log "串行: $serial_time | MPI-2: $mpi2_time | MPI-4: $mpi4_time" 