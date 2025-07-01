

CXX = g++
NVCC = nvcc
CXXFLAGS = -std=c++20 -O3 -Wall -Wextra
NVCCFLAGS = -std=c++20 -O3 -arch=sm_75

SRCDIR = src
INCDIR = src/include
TESTDIR = .

INCLUDES = -I$(INCDIR) -I$(INCDIR)/general -I$(INCDIR)/CUDA

all: main_cpu

main_cpu: main.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DUSE_CPU_ONLY -o $@ $<

main_gpu: main.cc $(INCDIR)/CUDA/ntt_gpu.o
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -o $@ $^

$(INCDIR)/CUDA/ntt_gpu.o: $(INCDIR)/CUDA/ntt_gpu.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -c -o $@ $<

perf: performance_compare.cu $(INCDIR)/CUDA/ntt_gpu.o
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -o performance_compare $^
	@echo "正在运行性能对比测试..."
	./performance_compare

clean:
	rm -f main_cpu main_gpu performance_compare $(INCDIR)/CUDA/ntt_gpu.o *.csv

gpu-info:
	@nvidia-smi || echo "NVIDIA GPU不可用"
	@nvcc --version || echo "NVCC编译器不可用"

help:
	@echo "Available targets:"
	@echo "  all       - 编译CPU版本 (默认)"
	@echo "  main_cpu  - 编译CPU版本"
	@echo "  main_gpu  - 编译GPU版本"
	@echo "  perf      - 编译并运行性能对比测试"
	@echo "  clean     - 清理编译产物"
	@echo "  gpu-info  - 显示GPU环境信息"
	@echo "  help      - 显示此帮助信息"

.PHONY: all clean gpu-info help perf 