cd ntt
g++ main.cc -o main -O2 -fopenmp -lpthread -std=c++11
if [ $? -eq 0 ]; then
  # 参数1：实验编号，pthread对应编号2，openmp对应编号3
  # 参数2：申请核心数
  # 参数3：申请线程数
  ./test.sh 2 1 1
else
  echo "编译失败，跳过执行 test.sh"
fi
