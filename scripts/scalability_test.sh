#!/bin/bash

# MPI NTT 可扩展性测试脚本
# 分析强可扩展性（固定问题规模，增加进程数）和弱可扩展性（按比例增加问题规模和进程数）

set -e

# 配置
RESULT_DIR="scalability_results"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${RESULT_DIR}/scalability_${TIMESTAMP}.csv"

# 创建目录
mkdir -p ${RESULT_DIR}

# CSV头部
echo "测试类型,进程数,问题规模,模数,执行时间(秒),加速比,效率" > ${LOG_FILE}

log() {
    echo "[$(date '+%H:%M:%S')] $1"
}

# 运行单次测试并返回时间（秒）
run_timed_test() {
    local np=$1
    local n=$2
    local mod=$3
    local description=$4
    
    log "测试: $description (np=$np, n=$n, mod=$mod)"
    
    cd ntt
    
    local start_time=$(date +%s.%N)
    
    echo "$n $mod" | mpirun -np $np ./main_mpi > /dev/null 2>&1
    
    local end_time=$(date +%s.%N)
    local duration=$(echo "$end_time - $start_time" | bc -l)
    
    cd ..
    
    echo $duration
}

# 强可扩展性测试
test_strong_scalability() {
    log "开始强可扩展性测试..."
    
    local n=131072
    local mod=104857601
    local process_counts=(1 2 4 8)
    local baseline_time=""
    
    for np in "${process_counts[@]}"; do
        local time_taken=$(run_timed_test $np $n $mod "强可扩展性-${np}进程")
        
        if [ $np -eq 1 ]; then
            baseline_time=$time_taken
            local speedup=1.0
            local efficiency=1.0
        else
            local speedup=$(echo "scale=3; $baseline_time / $time_taken" | bc -l)
            local efficiency=$(echo "scale=3; $speedup / $np" | bc -l)
        fi
        
        echo "强可扩展性,$np,$n,$mod,$time_taken,$speedup,$efficiency" >> ${LOG_FILE}
        
        log "结果: 时间=${time_taken}s, 加速比=${speedup}, 效率=${efficiency}"
    done
}

# 弱可扩展性测试（假设问题规模可以按进程数比例调整）
test_weak_scalability() {
    log "开始弱可扩展性测试..."
    
    local base_n=65536  # 基础问题规模
    local mod=104857601
    local process_counts=(1 2 4 8)
    local baseline_time=""
    
    for np in "${process_counts[@]}"; do
        local n=$((base_n * np))  # 按进程数比例增加问题规模
        local time_taken=$(run_timed_test $np $n $mod "弱可扩展性-${np}进程")
        
        if [ $np -eq 1 ]; then
            baseline_time=$time_taken
            local efficiency=1.0
        else
            local efficiency=$(echo "scale=3; $baseline_time / $time_taken" | bc -l)
        fi
        
        echo "弱可扩展性,$np,$n,$mod,$time_taken,N/A,$efficiency" >> ${LOG_FILE}
        
        log "结果: n=$n, 时间=${time_taken}s, 效率=${efficiency}"
    done
}

# 不同模数的可扩展性测试
test_different_moduli() {
    log "开始不同模数的可扩展性测试..."
    
    local n=131072
    local moduli=(7340033 104857601 469762049)
    local mod_names=("小模数" "中模数" "大模数")
    local process_counts=(1 2 4)
    
    for i in "${!moduli[@]}"; do
        local mod="${moduli[$i]}"
        local mod_name="${mod_names[$i]}"
        local baseline_time=""
        
        log "测试模数: $mod ($mod_name)"
        
        for np in "${process_counts[@]}"; do
            local time_taken=$(run_timed_test $np $n $mod "${mod_name}-${np}进程")
            
            if [ $np -eq 1 ]; then
                baseline_time=$time_taken
                local speedup=1.0
                local efficiency=1.0
            else
                local speedup=$(echo "scale=3; $baseline_time / $time_taken" | bc -l)
                local efficiency=$(echo "scale=3; $speedup / $np" | bc -l)
            fi
            
            echo "${mod_name}可扩展性,$np,$n,$mod,$time_taken,$speedup,$efficiency" >> ${LOG_FILE}
            
            log "结果: 时间=${time_taken}s, 加速比=${speedup}, 效率=${efficiency}"
        done
    done
}

# 混合并行可扩展性测试
test_hybrid_scalability() {
    log "开始混合并行可扩展性测试..."
    
    local n=131072
    local mod=104857601
    local configs=(
        "1 1"   # 1进程 1线程
        "2 1"   # 2进程 1线程
        "2 2"   # 2进程 2线程
        "4 1"   # 4进程 1线程
        "4 2"   # 4进程 2线程
    )
    local baseline_time=""
    
    for config in "${configs[@]}"; do
        local np=$(echo $config | cut -d' ' -f1)
        local threads=$(echo $config | cut -d' ' -f2)
        
        log "测试混合配置: ${np}进程 × ${threads}线程"
        
        cd ntt
        
        export OMP_NUM_THREADS=$threads
        local start_time=$(date +%s.%N)
        
        echo "$n $mod" | mpirun -np $np ./main_mpi > /dev/null 2>&1
        
        local end_time=$(date +%s.%N)
        local time_taken=$(echo "$end_time - $start_time" | bc -l)
        
        unset OMP_NUM_THREADS
        cd ..
        
        if [ "$config" = "1 1" ]; then
            baseline_time=$time_taken
            local speedup=1.0
            local efficiency=1.0
        else
            local total_cores=$((np * threads))
            local speedup=$(echo "scale=3; $baseline_time / $time_taken" | bc -l)
            local efficiency=$(echo "scale=3; $speedup / $total_cores" | bc -l)
        fi
        
        echo "混合并行,${np}×${threads},$n,$mod,$time_taken,$speedup,$efficiency" >> ${LOG_FILE}
        
        log "结果: 时间=${time_taken}s, 加速比=${speedup}, 效率=${efficiency}"
    done
}

# 生成可扩展性报告
generate_scalability_report() {
    log "生成可扩展性分析报告..."
    
    local report_file="${RESULT_DIR}/scalability_report_${TIMESTAMP}.md"
    
    cat > ${report_file} << 'EOF'
# MPI NTT 可扩展性分析报告

## 测试配置
- 测试时间: $(date)
- 问题规模: 主要测试 n=131072
- 主要模数: 104857601
- 进程数范围: 1-8

## 强可扩展性分析

强可扩展性测试固定问题规模，增加进程数，观察加速比和效率变化。

理想情况下：
- 加速比 = 进程数
- 效率 = 1.0 (100%)

## 弱可扩展性分析

弱可扩展性测试按比例增加问题规模和进程数，观察执行时间是否保持稳定。

理想情况下：
- 执行时间保持不变
- 效率 = 1.0 (100%)

## 混合并行分析

测试MPI进程数与OpenMP线程数的不同组合，寻找最优配置。

## 主要发现

1. **最优配置**：[待填充]
2. **可扩展性瓶颈**：[待填充]
3. **通信开销影响**：[待填充]

## 数据文件

详细数据请参考: scalability_TIMESTAMP.csv
EOF

    # 替换时间戳
    sed -i "s/TIMESTAMP/${TIMESTAMP}/g" ${report_file}
    
    log "报告已生成: ${report_file}"
}

# 主函数
main() {
    log "开始可扩展性测试..."
    log "结果将保存到: ${RESULT_DIR}"
    
    # 检查环境
    if [ ! -d "ntt" ]; then
        echo "错误: 请在项目根目录运行此脚本"
        exit 1
    fi
    
    if ! command -v mpirun &> /dev/null; then
        echo "错误: MPI环境未找到"
        exit 1
    fi
    
    if ! command -v bc &> /dev/null; then
        echo "警告: bc命令未找到，计算可能不准确"
    fi
    
    # 编译代码
    cd ntt
    log "编译代码..."
    mpic++ -O3 -fopenmp -DUSE_MPI -o main_mpi main.cc
    cd ..
    
    # 运行测试
    test_strong_scalability
    test_weak_scalability
    test_different_moduli
    test_hybrid_scalability
    
    # 生成报告
    generate_scalability_report
    
    log "可扩展性测试完成！"
    log "查看结果: ls -la ${RESULT_DIR}/"
}

# 运行主函数
main "$@" 