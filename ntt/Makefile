# GPU NTT并行化项目主Makefile
# 支持编译不同版本的GPU NTT实现

# 编译器设置
NVCC = nvcc
CXX = g++
CUDA_ARCH = -arch=sm_89
CUDA_FLAGS = -O3 -std=c++17 $(CUDA_ARCH)
CXX_FLAGS = -O3 -std=c++17

# 目录设置
SRC_DIR = src
CUDA_DIR = $(SRC_DIR)/cuda
INCLUDE_DIR = $(SRC_DIR)/include
TEST_DIR = $(SRC_DIR)/test
BUILD_DIR = build

# 源文件
CUDA_COMMON = $(CUDA_DIR)/gpu_ntt_common.h $(CUDA_DIR)/memory_wrapper.h
CUDA_V1 = $(CUDA_DIR)/gpu_ntt_v1.cu
CUDA_V2 = $(CUDA_DIR)/gpu_ntt_v2.cu
CUDA_V3 = $(CUDA_DIR)/gpu_ntt_v3.cu

# 目标文件
MAIN_TARGET = $(BUILD_DIR)/main
PERF_TARGET = $(BUILD_DIR)/simple_test

# 默认目标：编译v2版本（推荐）
all: $(MAIN_TARGET)

# 编译主程序（使用v2版本）
$(MAIN_TARGET): main.cc $(CUDA_V2) $(CUDA_COMMON)
	@mkdir -p $(BUILD_DIR)
	$(NVCC) $(CUDA_FLAGS) -I$(INCLUDE_DIR) -I$(CUDA_DIR) -o $@ main.cc $(CUDA_V2)

# 编译v1版本
v1: $(BUILD_DIR)/main_v1
$(BUILD_DIR)/main_v1: main.cc $(CUDA_V1) $(CUDA_COMMON)
	@mkdir -p $(BUILD_DIR)
	$(NVCC) $(CUDA_FLAGS) -DGPU_NTT_V1 -I$(INCLUDE_DIR) -I$(CUDA_DIR) -o $@ main.cc $(CUDA_V1)

# 编译v2修复版本
v2-fixed: $(BUILD_DIR)/main_v2_fixed
$(BUILD_DIR)/main_v2_fixed: main.cc $(CUDA_DIR)/gpu_ntt_v2_fixed.cu $(CUDA_COMMON)
	@mkdir -p $(BUILD_DIR)
	$(NVCC) $(CUDA_FLAGS) -DGPU_NTT_V2_FIXED -I$(INCLUDE_DIR) -I$(CUDA_DIR) -o $@ main.cc $(CUDA_DIR)/gpu_ntt_v2_fixed.cu

# 编译v2增强版本
v2-enhanced: $(BUILD_DIR)/main_v2_enhanced
$(BUILD_DIR)/main_v2_enhanced: main.cc $(CUDA_DIR)/gpu_ntt_v2_enhanced.cu $(CUDA_COMMON)
	@mkdir -p $(BUILD_DIR)
	$(NVCC) $(CUDA_FLAGS) -DGPU_NTT_V2_ENHANCED -I$(INCLUDE_DIR) -I$(CUDA_DIR) -o $@ main.cc $(CUDA_DIR)/gpu_ntt_v2_enhanced.cu

# 编译v3版本
v3: $(BUILD_DIR)/main_v3
$(BUILD_DIR)/main_v3: main.cc $(CUDA_V3) $(CUDA_COMMON)
	@mkdir -p $(BUILD_DIR)
	$(NVCC) $(CUDA_FLAGS) -DGPU_NTT_V3 -I$(INCLUDE_DIR) -I$(CUDA_DIR) -o $@ main.cc $(CUDA_V3)

# 编译性能测试
test: $(PERF_TARGET)
$(PERF_TARGET): $(TEST_DIR)/simple_test.cc $(CUDA_V2) $(CUDA_COMMON)
	@mkdir -p $(BUILD_DIR)
	$(NVCC) $(CUDA_FLAGS) -I$(INCLUDE_DIR) -I$(CUDA_DIR) -o $@ $(TEST_DIR)/simple_test.cc $(CUDA_V2)

# 运行测试
run-test: $(PERF_TARGET)
	./$(PERF_TARGET)

# 清理编译产物
clean:
	rm -rf $(BUILD_DIR)
	find . -name "*.o" -delete
	find . -name "*.tmp" -delete

# 深度清理（包括旧文件）
distclean: clean
	rm -f main_cuda final_performance_benchmark test_v3_*
	rm -f *.o *.tmp

# 显示帮助信息
help:
	@echo "GPU NTT并行化项目编译选项："
	@echo "  make        - 编译主程序（v2版本，推荐）"
	@echo "  make v1     - 编译v1版本（预计算旋转因子优化）"
	@echo "  make v2-fixed - 编译v2修复版本（安全共享内存优化）"
	@echo "  make v2-enhanced - 编译v2增强版本（性能优化）"
	@echo "  make v3     - 编译v3版本（多层共享内存实验）"
	@echo "  make test   - 编译性能测试程序"
	@echo "  make run-test - 运行性能测试"
	@echo "  make clean  - 清理编译产物"
	@echo "  make distclean - 深度清理所有文件"
	@echo "  make help   - 显示此帮助信息"

.PHONY: all v1 v3 test run-test clean distclean help 