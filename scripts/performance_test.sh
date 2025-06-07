#!/bin/bash

# MPI NTT性能测试脚本
# 作者：Lab4实验组
# 功能：测试不同问题规模、进程数、混合并行配置下的性能

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 测试配置
RESULT_DIR="performance_results"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${RESULT_DIR}/performance_test_${TIMESTAMP}.log"

# 创建结果目录
mkdir -p ${RESULT_DIR}

# 日志函数
log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] $1${NC}" | tee -a ${LOG_FILE}
}

warn() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $1${NC}" | tee -a ${LOG_FILE}
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1${NC}" | tee -a ${LOG_FILE}
}

# 检查环境
check_environment() {
    log "检查测试环境..."
    
    # 检查MPI环境
    if ! command -v mpirun &> /dev/null; then
        error "MPI环境未找到，请安装OpenMPI或MPICH"
        exit 1
    fi
    
    # 检查编译器
    if ! command -v mpic++ &> /dev/null; then
        error "MPI C++编译器未找到"
        exit 1
    fi
    
    # 检查代码目录
    if [ ! -d "ntt" ]; then
        error "ntt代码目录未找到"
        exit 1
    fi
    
    log "环境检查通过"
}

# 编译代码
compile_code() {
    log "编译测试代码..."
    
    cd ntt
    
    # 注意：由于代码结构原因，我们统一使用MPI版本，通过进程数控制并行度
    log "编译代码..."
    
    # 编译MPI版本
    log "编译MPI版本..."
    mpic++ -O3 -fopenmp -DUSE_MPI -o main_mpi main.cc
    
    cd ..
    log "编译完成"
}

# 单次测试函数
run_single_test() {
    local test_name="$1"
    local command="$2"
    local description="$3"
    
    echo -e "${BLUE}=== $test_name ===${NC}" | tee -a ${LOG_FILE}
    echo "描述: $description" | tee -a ${LOG_FILE}
    echo "命令: $command" | tee -a ${LOG_FILE}
    echo "开始时间: $(date)" | tee -a ${LOG_FILE}
    
    # 运行测试并记录时间
    local start_time=$(date +%s.%N)
    
    # 设置超时时间（5分钟）
    timeout 300 bash -c "$command" >> ${LOG_FILE} 2>&1
    local exit_code=$?
    
    local end_time=$(date +%s.%N)
    local duration=$(echo "$end_time - $start_time" | bc -l)
    
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ 测试完成，耗时: ${duration}秒${NC}" | tee -a ${LOG_FILE}
    elif [ $exit_code -eq 124 ]; then
        echo -e "${YELLOW}⚠ 测试超时（5分钟）${NC}" | tee -a ${LOG_FILE}
    else
        echo -e "${RED}✗ 测试失败，退出码: $exit_code${NC}" | tee -a ${LOG_FILE}
    fi
    
    echo "" | tee -a ${LOG_FILE}
    return $exit_code
}

# 测试1：串行性能基准测试
test_serial_baseline() {
    log "开始串行性能基准测试..."
    
    cd ntt
    
    # 小规模测试（验证正确性）
    run_single_test "串行-小规模" \
        "./main_serial" \
        "n=4的小规模测试，验证算法正确性"
    
    # 大规模测试（性能基准）
    run_single_test "串行-大规模-模数1" \
        "./main_serial < ../data/input1.txt" \
        "n=131072，模数=7340033的串行性能基准"
    
    run_single_test "串行-大规模-模数2" \
        "./main_serial < ../data/input2.txt" \
        "n=131072，模数=104857601的串行性能基准"
    
    run_single_test "串行-大规模-模数3" \
        "./main_serial < ../data/input3.txt" \
        "n=131072，模数=469762049的串行性能基准"
    
    cd ..
}

# 测试2：MPI可扩展性测试
test_mpi_scalability() {
    log "开始MPI可扩展性测试..."
    
    cd ntt
    
    local process_counts=(1 2 4 8)
    local test_cases=("input1.txt" "input2.txt" "input3.txt")
    local descriptions=("模数=7340033" "模数=104857601" "模数=469762049")
    
    for i in "${!test_cases[@]}"; do
        local input_file="${test_cases[$i]}"
        local desc="${descriptions[$i]}"
        
        log "测试输入文件: $input_file ($desc)"
        
        for np in "${process_counts[@]}"; do
            run_single_test "MPI-${np}进程-${desc}" \
                "mpirun -np $np ./main_mpi < ../data/$input_file" \
                "使用${np}个进程测试$desc"
        done
    done
    
    cd ..
}

# 测试3：混合并行测试（MPI + OpenMP）
test_hybrid_parallel() {
    log "开始混合并行测试..."
    
    cd ntt
    
    # 设置不同的OpenMP线程数
    local omp_threads=(1 2 4)
    local mpi_processes=(2 4)
    
    for np in "${mpi_processes[@]}"; do
        for threads in "${omp_threads[@]}"; do
            export OMP_NUM_THREADS=$threads
            
            run_single_test "混合-${np}进程${threads}线程" \
                "mpirun -np $np ./main_mpi < ../data/input2.txt" \
                "MPI ${np}进程 + OpenMP ${threads}线程"
        done
    done
    
    # 重置OpenMP线程数
    unset OMP_NUM_THREADS
    
    cd ..
}

# 测试4：通信开销分析
test_communication_overhead() {
    log "开始通信开销分析测试..."
    
    cd ntt
    
    # 使用调试版本来获取通信开销信息
    log "编译调试版本..."
    mpic++ -O3 -fopenmp -DUSE_MPI -DDEBUG_TIMING -o main_mpi_debug main.cc
    
    # 测试不同进程数下的通信开销
    local process_counts=(2 4 8)
    
    for np in "${process_counts[@]}"; do
        run_single_test "通信开销-${np}进程" \
            "mpirun -np $np ./main_mpi_debug < ../data/input2.txt" \
            "分析${np}进程下的通信开销"
    done
    
    cd ..
}

# 测试5：不同问题规模的强可扩展性测试
test_strong_scalability() {
    log "开始强可扩展性测试..."
    
    cd ntt
    
    # 固定问题规模，增加进程数
    local process_counts=(1 2 4 8 16)
    
    for np in "${process_counts[@]}"; do
        if [ $np -eq 1 ]; then
            run_single_test "强可扩展性-${np}进程" \
                "./main_serial < ../data/input2.txt" \
                "强可扩展性基准：单进程串行版本"
        else
            run_single_test "强可扩展性-${np}进程" \
                "mpirun -np $np ./main_mpi < ../data/input2.txt" \
                "强可扩展性测试：${np}进程"
        fi
    done
    
    cd ..
}

# 测试6：性能profile分析
test_performance_profiling() {
    log "开始性能profile分析..."
    
    cd ntt
    
    # 使用time命令进行详细计时
    run_single_test "Profile-串行" \
        "/usr/bin/time -v ./main_serial < ../data/input2.txt" \
        "串行版本详细性能profile"
    
    run_single_test "Profile-MPI-2进程" \
        "/usr/bin/time -v mpirun -np 2 ./main_mpi < ../data/input2.txt" \
        "MPI 2进程详细性能profile"
    
    run_single_test "Profile-MPI-4进程" \
        "/usr/bin/time -v mpirun -np 4 ./main_mpi < ../data/input2.txt" \
        "MPI 4进程详细性能profile"
    
    cd ..
}

# 生成性能报告
generate_performance_report() {
    log "生成性能测试报告..."
    
    local report_file="${RESULT_DIR}/performance_report_${TIMESTAMP}.md"
    
    cat > ${report_file} << EOF
# MPI NTT 性能测试报告

## 测试环境
- 测试时间: $(date)
- 系统信息: $(uname -a)
- 处理器信息: $(cat /proc/cpuinfo | grep "model name" | head -1 | cut -d: -f2)
- 内存信息: $(free -h | grep Mem)
- MPI版本: $(mpirun --version | head -1)

## 测试配置
- 编译器: $(mpic++ --version | head -1)
- 编译选项: -O3 -fopenmp -DUSE_MPI
- 测试问题规模: n=131072
- 测试模数: 7340033, 104857601, 469762049

## 测试结果

详细的测试日志请参考: ${LOG_FILE}

### 主要性能指标

1. **串行基准性能**
   - 待填充具体数据

2. **MPI可扩展性**
   - 待填充具体数据

3. **混合并行效果**
   - 待填充具体数据

4. **通信开销分析**
   - 待填充具体数据

## 结论

基于测试结果，我们可以得出以下结论：
1. [待补充]
2. [待补充]
3. [待补充]

EOF

    log "性能报告已生成: ${report_file}"
}

# 主函数
main() {
    log "开始MPI NTT性能测试"
    log "测试结果将保存到: ${RESULT_DIR}"
    log "详细日志: ${LOG_FILE}"
    
    # 检查环境
    check_environment
    
    # 编译代码
    compile_code
    
    # 运行各项测试
    test_serial_baseline
    test_mpi_scalability
    test_hybrid_parallel
    test_communication_overhead
    test_strong_scalability
    test_performance_profiling
    
    # 生成报告
    generate_performance_report
    
    log "所有性能测试完成！"
    log "查看结果: ls -la ${RESULT_DIR}/"
}

# 检查bc命令（用于浮点运算）
if ! command -v bc &> /dev/null; then
    warn "bc命令未找到，时间计算可能不准确"
fi

# 运行主函数
main "$@" 