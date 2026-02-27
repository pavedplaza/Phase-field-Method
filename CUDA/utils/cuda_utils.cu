/**
 * @file cuda_utils.cu
 * @brief CUDA工具函数实现
 * @author pavedplaza
 * @date 2025-02-15
 *
 * 实现GPU信息查询、内存管理等辅助函数
 */

#include "cuda_headers.h"
#include "phase_field_params.h"

/**
 * @brief 打印GPU设备信息
 *
 * 显示所有可用GPU的详细信息，包括：
 * - 设备名称和计算能力
 * - 全局内存大小
 * - Shared memory大小
 * - 最大线程数等
 */
void print_gpu_info() {
    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));

    if (device_count == 0) {
        fprintf(stderr, "错误：未检测到CUDA设备！\n");
        exit(EXIT_FAILURE);
    }

    printf("==============================================\n");
    printf("检测到 %d 个CUDA设备\n", device_count);
    printf("==============================================\n");

    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, i));

        printf("\n[GPU %d] %s\n", i, prop.name);
        printf("  -------------------------------------------\n");
        printf("  计算能力:           %d.%d\n", prop.major, prop.minor);
        printf("  全局内存:           %.2f GB\n", prop.totalGlobalMem / 1e9);
        printf("  Shared Memory/Block: %zu KB\n", prop.sharedMemPerBlock / 1024);
        printf("  最大线程/Block:     %d\n", prop.maxThreadsPerBlock);
        printf("  最大线程维度:       [%d, %d, %d]\n",
               prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
        printf("  最大Block维度:      [%d, %d, %d]\n",
               prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
        printf("  Warp大小:           %d\n", prop.warpSize);
        printf("  SM数量:             %d\n", prop.multiProcessorCount);
        printf("  L2缓存大小:         %d KB\n", prop.l2CacheSize / 1024);
        // clockRate 在 CUDA 11.0+ 中已被移除，不再显示
        // printf("  最大时钟频率:       %.1f GHz\n", prop.clockRate / 1e6);

        // 计算理论性能（使用时钟频率估算值，或者跳过）
        // float tflops = (prop.multiProcessorCount * prop.clockRate * 1e-6f *
        //                prop.major * 128.f) / 1e12f * 2.f;
        // printf("  理论峰值性能:       ~%.1f TFLOPS\n", tflops);
        printf("  理论峰值性能:       (需要查询GPU规格)\n");
    }

    printf("\n==============================================\n");
    printf("当前使用设备: ");
    int current_device = 0;
    CUDA_CHECK(cudaGetDevice(&current_device));
    cudaDeviceProp current_prop;
    CUDA_CHECK(cudaGetDeviceProperties(&current_prop, current_device));
    printf("%s (GPU %d)\n", current_prop.name, current_device);
    printf("==============================================\n\n");
}

/**
 * @brief 初始化GPU数组
 *
 * 在GPU上分配所有必要的数组内存
 *
 * @param arrays GPU数组结构体
 * @param params 模拟参数
 * @return 成功返回true，失败返回false
 */
bool init_gpu_arrays(GPUArrays& arrays, const PhaseFieldParams& params) {
    printf("初始化GPU数组...\n");

    size_t total_size = params.Nx * params.Ny * sizeof(float);
    size_t total_allocated = 0;

    // 分配相场数组
    CUDA_CHECK(cudaMalloc(&arrays.d_phi, total_size));
    CUDA_CHECK(cudaMalloc(&arrays.d_phi_new, total_size));
    total_allocated += 2 * total_size;

    // 分配浓度场数组
    CUDA_CHECK(cudaMalloc(&arrays.d_U, total_size));
    CUDA_CHECK(cudaMalloc(&arrays.d_U_new, total_size));
    total_allocated += 2 * total_size;

    // Phase 2: 分配梯度数组
    CUDA_CHECK(cudaMalloc(&arrays.d_phidx, total_size));
    CUDA_CHECK(cudaMalloc(&arrays.d_phidy, total_size));
    total_allocated += 2 * total_size;

    // Phase 2: 分配各向异性参数数组
    CUDA_CHECK(cudaMalloc(&arrays.d_epsilon, total_size));
    CUDA_CHECK(cudaMalloc(&arrays.d_epsilon_deriv, total_size));
    total_allocated += 2 * total_size;

    // Phase 3: 分配浓度场梯度数组
    CUDA_CHECK(cudaMalloc(&arrays.d_Udx, total_size));
    CUDA_CHECK(cudaMalloc(&arrays.d_Udy, total_size));
    CUDA_CHECK(cudaMalloc(&arrays.d_phi_change_rate, total_size));
    total_allocated += 3 * total_size;

    // 保存数组尺寸信息
    arrays.Nx = params.Nx;
    arrays.Ny = params.Ny;

    printf("  GPU内存分配完成: %.2f MB\n", total_allocated / (1024.0 * 1024.0));

    return true;
}

/**
 * @brief 释放GPU数组
 *
 * 释放所有GPU内存
 *
 * @param arrays GPU数组结构体
 */
void free_gpu_arrays(GPUArrays& arrays) {
    printf("释放GPU数组...\n");

    if (arrays.d_phi != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_phi));
        arrays.d_phi = nullptr;
    }
    if (arrays.d_phi_new != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_phi_new));
        arrays.d_phi_new = nullptr;
    }
    if (arrays.d_U != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_U));
        arrays.d_U = nullptr;
    }
    if (arrays.d_U_new != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_U_new));
        arrays.d_U_new = nullptr;
    }
    if (arrays.d_phidx != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_phidx));
        arrays.d_phidx = nullptr;
    }
    if (arrays.d_phidy != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_phidy));
        arrays.d_phidy = nullptr;
    }
    if (arrays.d_epsilon != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_epsilon));
        arrays.d_epsilon = nullptr;
    }
    if (arrays.d_epsilon_deriv != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_epsilon_deriv));
        arrays.d_epsilon_deriv = nullptr;
    }
    // Phase 3: 释放浓度场梯度数组
    if (arrays.d_Udx != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_Udx));
        arrays.d_Udx = nullptr;
    }
    if (arrays.d_Udy != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_Udy));
        arrays.d_Udy = nullptr;
    }
    if (arrays.d_phi_change_rate != nullptr) {
        CUDA_CHECK(cudaFree(arrays.d_phi_change_rate));
        arrays.d_phi_change_rate = nullptr;
    }

    printf("  GPU内存释放完成\n");
}

/**
 * @brief 检查GPU内存是否足够
 *
 * 在分配内存前检查是否有足够的GPU内存
 *
 * @param required_bytes 需要的字节数
 * @return 内存足够返回true，否则返回false
 */
bool check_gpu_memory(size_t required_bytes) {
    size_t free_bytes = 0, total_bytes = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_bytes, &total_bytes));

    printf("GPU内存状态: %.2f GB 可用 / %.2f GB 总计\n",
           free_bytes / 1e9, total_bytes / 1e9);

    if (required_bytes > free_bytes) {
        fprintf(stderr, "错误：需要 %.2f GB 内存，但只有 %.2f GB 可用\n",
                required_bytes / 1e9, free_bytes / 1e9);
        return false;
    }

    return true;
}

/**
 * @brief 重置GPU设备
 *
 * 清理GPU状态，释放所有资源
 * 主要用于调试和错误恢复
 */
void reset_gpu_device() {
    printf("重置GPU设备...\n");
    CUDA_CHECK(cudaDeviceReset());
}

/**
 * @brief 设置GPU设备
 *
 * 选择要使用的GPU设备
 *
 * @param device_id 设备ID（从0开始）
 */
void set_gpu_device(int device_id) {
    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));

    if (device_id < 0 || device_id >= device_count) {
        fprintf(stderr, "错误：设备ID %d 无效（范围: 0-%d）\n",
                device_id, device_count - 1);
        exit(EXIT_FAILURE);
    }

    CUDA_CHECK(cudaSetDevice(device_id));
    printf("使用GPU设备 %d\n", device_id);
}
