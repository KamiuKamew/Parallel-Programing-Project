cd ntt
g++ main.cc -o main -O2 -fopenmp -lpthread -std=c++11
if [ $? -eq 0 ]; then
  ./main
else
  echo "编译失败，跳过执行 test.sh"
fi
