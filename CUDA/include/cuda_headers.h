/**
 * @file cuda_headers.h
 * @brief CUDA通用头文件和工具宏
 * @author PhaseFieldCUDA Team
 * @date 2025-02-15
 *
 * 包含所有CUDA相关的头文件和常用宏定义
 */

#pragma once

// ========== CUDA头文件 ==========
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// ========== C++标准库 ==========
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

// ========== CUDA错误检查宏 ==========

/**
 * @brief CUDA错误检查宏
 *
 * 用法：CUDA_CHECK(cudaMalloc(&ptr, size));
 *
 * 如果CUDA调用失败，打印错误信息并退出程序
 */
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = (call); \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s (error code: %d)\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err), err); \
            fflush(stderr); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

/**
 * @brief CUDA kernel启动错误检查宏
 *
 * 用法：CUDA_KERNEL_CHECK();
 *
 * 检查kernel启动是否有错误（异步错误）
 */
#define CUDA_KERNEL_CHECK() \
    do { \
        cudaError_t err = cudaGetLastError(); \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA kernel error at %s:%d: %s (error code: %d)\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err), err); \
            fflush(stderr); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

/**
 * @brief CUDA设备同步和错误检查宏
 *
 * 用法：CUDA_SYNC_CHECK();
 *
 * 同步设备并检查是否有错误
 */
#define CUDA_SYNC_CHECK() \
    do { \
        CUDA_CHECK(cudaDeviceSynchronize()); \
        CUDA_KERNEL_CHECK(); \
    } while(0)

// ========== CUDA数学常量 ==========

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_F
#define M_PI_F 3.14159265358979323846f
#endif

// ========== 常用CUDA设备函数 ==========

/**
 * @brief 设备函数：计算周期性边界索引
 *
 * @param i 当前索引
 * @param max 最大索引值
 * @return 周期性边界索引
 */
__device__ __forceinline__ int get_periodic_index(int i, int max) {
    return (i + max) % max;
}

/**
 * @brief 设备函数：安全的浮点除法（避免除零）
 *
 * @param numerator 分子
 * @param denominator 分母
 * @param default_value 分母为零时的默认返回值
 * @return 除法结果或默认值
 */
__device__ __forceinline__ float safe_divide(float numerator,
                                              float denominator,
                                              float default_value = 0.0f) {
    return (fabsf(denominator) > 1e-12f) ? (numerator / denominator) : default_value;
}

/**
 * @brief 设备函数：计算梯度幅值
 *
 * @param gx x方向梯度
 * @param gy y方向梯度
 * @return 梯度幅值
 */
__device__ __forceinline__ float gradient_magnitude(float gx, float gy) {
    return sqrtf(gx * gx + gy * gy);
}

/**
 * @brief 设备函数：计算角度（atan2的包装）
 *
 * @param y y坐标
 * @param x x坐标
 * @return 角度（弧度），范围[-π, π]
 */
__device__ __forceinline__ float compute_angle(float y, float x) {
    return atan2f(y, x);
}

// ========== 性能计时工具 ==========

#ifdef __CUDACC__

/**
 * @brief CUDA事件计时器类
 *
 * 用于精确测量GPU kernel执行时间
 */
class CudaTimer {
public:
    CudaTimer() {
        CUDA_CHECK(cudaEventCreate(&start_));
        CUDA_CHECK(cudaEventCreate(&stop_));
    }

    ~CudaTimer() {
        CUDA_CHECK(cudaEventDestroy(start_));
        CUDA_CHECK(cudaEventDestroy(stop_));
    }

    void start() {
        CUDA_CHECK(cudaEventRecord(start_));
    }

    void stop() {
        CUDA_CHECK(cudaEventRecord(stop_));
        CUDA_CHECK(cudaEventSynchronize(stop_));
    }

    float elapsed_milliseconds() {
        float ms = 0.0f;
        CUDA_CHECK(cudaEventElapsedTime(&ms, start_, stop_));
        return ms;
    }

private:
    cudaEvent_t start_;
    cudaEvent_t stop_;
};

#endif // __CUDACC__

// ========== 内存管理辅助函数 ==========

/**
 * @brief 安全的GPU内存分配
 *
 * @param ptr 指针指针
 * @param size 字节数
 * @param name 分配的内存名称（用于错误信息）
 */
inline void cuda_malloc_safe(void** ptr, size_t size, const char* name = "unknown") {
    CUDA_CHECK(cudaMalloc(ptr, size));
    printf("Allocated %.2f MB for %s\n", size / (1024.0 * 1024.0), name);
}

/**
 * @brief 安全的GPU内存释放
 *
 * @param ptr 要释放的指针
 * @param name 释放的内存名称（用于错误信息）
 */
inline void cuda_free_safe(void* ptr, const char* name = "unknown") {
    if (ptr != nullptr) {
        CUDA_CHECK(cudaFree(ptr));
        printf("Freed memory for %s\n", name);
    }
}
