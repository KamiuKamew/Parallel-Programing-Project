#include "src/include/CUDA/ntt.h"
#include "src/include/ntt.h"
#include "src/include/transform.h"
#include <cstring>
#include <iostream>
#include <vector>

void print_array(const char *name, u64 *arr, u64 n) {
  printf("%s: ", name);
  for (u64 i = 0; i < n; i++) {
    printf("%lu ", arr[i]);
  }
  printf("\n");
}

int main() {
  printf("=== 小规模GPU NTT调试测试 ===\n\n");

  // 使用小规模数据: n=8
  u64 n = 8;
  u64 p = 7340033;
  u64 omega = 3;

  // 简单测试数据
  u64 a_cpu[8] = {1, 2, 3, 4, 0, 0, 0, 0};
  u64 a_gpu[8];
  memcpy(a_gpu, a_cpu, sizeof(a_cpu));

  printf("原始数据:\n");
  print_array("数据", a_cpu, n);

  printf("\n--- CPU NTT正变换 ---\n");
  ntt_forward(a_cpu, n, p, omega);
  print_array("CPU正变换结果", a_cpu, n);

  printf("\n--- GPU NTT正变换 ---\n");
  ntt_forward_gpu(a_gpu, n, p, omega);
  print_array("GPU正变换结果", a_gpu, n);

  // 验证正变换结果
  printf("\n--- 正变换结果对比 ---\n");
  bool forward_match = true;
  for (u64 i = 0; i < n; i++) {
    if (a_cpu[i] != a_gpu[i]) {
      printf("位置%lu不匹配: CPU=%lu, GPU=%lu\n", i, a_cpu[i], a_gpu[i]);
      forward_match = false;
    }
  }
  if (forward_match) {
    printf("✓ 正变换结果一致\n");
  } else {
    printf("✗ 正变换结果不一致\n");
  }

  // 测试逆变换
  printf("\n--- CPU NTT逆变换 ---\n");
  ntt_inverse(a_cpu, n, p, omega);
  print_array("CPU逆变换结果", a_cpu, n);

  printf("\n--- GPU NTT逆变换 ---\n");
  ntt_inverse_gpu(a_gpu, n, p, omega);
  print_array("GPU逆变换结果", a_gpu, n);

  // 验证逆变换结果
  printf("\n--- 逆变换结果对比 ---\n");
  bool inverse_match = true;
  for (u64 i = 0; i < n; i++) {
    if (a_cpu[i] != a_gpu[i]) {
      printf("位置%lu不匹配: CPU=%lu, GPU=%lu\n", i, a_cpu[i], a_gpu[i]);
      inverse_match = false;
    }
  }
  if (inverse_match) {
    printf("✓ 逆变换结果一致\n");
  } else {
    printf("✗ 逆变换结果不一致\n");
  }

  return 0;
}