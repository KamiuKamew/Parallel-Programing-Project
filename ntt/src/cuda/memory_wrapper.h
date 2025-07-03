#ifndef MEMORY_SAFE_WRAPPER_H
#define MEMORY_SAFE_WRAPPER_H

#include <cuda_runtime.h>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cstdint>

// 类型定义
using u64 = uint64_t;

// GPU内存的RAII包装器
template <typename T>
class CudaMemoryWrapper
{
private:
    T *ptr_;
    size_t size_;

public:
    // 构造函数：分配GPU内存
    explicit CudaMemoryWrapper(size_t count) : size_(count)
    {
        cudaError_t err = cudaMalloc((void **)&ptr_, count * sizeof(T));
        if (err != cudaSuccess)
        {
            throw std::runtime_error("GPU内存分配失败: " + std::string(cudaGetErrorString(err)));
        }
        std::cout << "✅ GPU内存分配成功: " << count * sizeof(T) << " bytes" << std::endl;
    }

    // 禁用拷贝构造和赋值
    CudaMemoryWrapper(const CudaMemoryWrapper &) = delete;
    CudaMemoryWrapper &operator=(const CudaMemoryWrapper &) = delete;

    // 移动构造函数
    CudaMemoryWrapper(CudaMemoryWrapper &&other) noexcept : ptr_(other.ptr_), size_(other.size_)
    {
        other.ptr_ = nullptr;
        other.size_ = 0;
    }

    // 移动赋值操作符
    CudaMemoryWrapper &operator=(CudaMemoryWrapper &&other) noexcept
    {
        if (this != &other)
        {
            cleanup();
            ptr_ = other.ptr_;
            size_ = other.size_;
            other.ptr_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }

    // 析构函数：自动释放GPU内存
    ~CudaMemoryWrapper()
    {
        cleanup();
    }

    // 获取原始指针
    T *get() const { return ptr_; }

    // 获取大小
    size_t size() const { return size_; }

    // 从CPU复制数据到GPU
    void copyFromHost(const T *host_data, size_t count)
    {
        if (count > size_)
        {
            throw std::runtime_error("复制数据超过分配的GPU内存大小");
        }
        cudaError_t err = cudaMemcpy(ptr_, host_data, count * sizeof(T), cudaMemcpyHostToDevice);
        if (err != cudaSuccess)
        {
            throw std::runtime_error("CPU到GPU数据复制失败: " + std::string(cudaGetErrorString(err)));
        }
    }

    // 从GPU复制数据到CPU
    void copyToHost(T *host_data, size_t count) const
    {
        if (count > size_)
        {
            throw std::runtime_error("复制数据超过分配的GPU内存大小");
        }
        cudaError_t err = cudaMemcpy(host_data, ptr_, count * sizeof(T), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess)
        {
            throw std::runtime_error("GPU到CPU数据复制失败: " + std::string(cudaGetErrorString(err)));
        }
    }

private:
    void cleanup()
    {
        if (ptr_)
        {
            cudaFree(ptr_);
            std::cout << "✅ GPU内存释放成功" << std::endl;
            ptr_ = nullptr;
            size_ = 0;
        }
    }
};

// CPU内存的RAII包装器
template <typename T>
class CpuMemoryWrapper
{
private:
    std::unique_ptr<T[]> ptr_;
    size_t size_;

public:
    // 构造函数：分配CPU内存
    explicit CpuMemoryWrapper(size_t count) : size_(count)
    {
        ptr_ = std::make_unique<T[]>(count);
        std::cout << "✅ CPU内存分配成功: " << count * sizeof(T) << " bytes" << std::endl;
    }

    // 析构函数自动处理（unique_ptr）
    ~CpuMemoryWrapper()
    {
        std::cout << "✅ CPU内存自动释放" << std::endl;
    }

    // 禁用拷贝构造和赋值
    CpuMemoryWrapper(const CpuMemoryWrapper &) = delete;
    CpuMemoryWrapper &operator=(const CpuMemoryWrapper &) = delete;

    // 允许移动
    CpuMemoryWrapper(CpuMemoryWrapper &&) = default;
    CpuMemoryWrapper &operator=(CpuMemoryWrapper &&) = default;

    // 获取原始指针
    T *get() const { return ptr_.get(); }

    // 获取大小
    size_t size() const { return size_; }

    // 数组访问操作符
    T &operator[](size_t index) { return ptr_[index]; }
    const T &operator[](size_t index) const { return ptr_[index]; }

    // 初始化所有元素为零
    void zero()
    {
        std::fill(ptr_.get(), ptr_.get() + size_, T{0});
    }

    // 复制数据从另一个数组
    void copyFrom(const T *source, size_t count)
    {
        if (count > size_)
        {
            throw std::runtime_error("复制数据超过分配的CPU内存大小");
        }
        std::copy(source, source + count, ptr_.get());
    }
};

// 便利的类型定义
using CudaU64Memory = CudaMemoryWrapper<u64>;
using CpuU64Memory = CpuMemoryWrapper<u64>;

// 辅助函数：检查CUDA错误
inline void checkCudaError(const std::string &operation)
{
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        throw std::runtime_error(operation + " 失败: " + std::string(cudaGetErrorString(err)));
    }
}

// 辅助函数：同步GPU并检查错误
inline void syncAndCheck(const std::string &operation)
{
    cudaDeviceSynchronize();
    checkCudaError(operation);
}

#endif // MEMORY_SAFE_WRAPPER_H