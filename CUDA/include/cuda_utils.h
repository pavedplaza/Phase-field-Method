/**
 * @file cuda_utils.h
 * @brief CUDA工具函数声明
 * @author PhaseFieldCUDA Team
 * @date 2025-02-26
 *
 * 声明CUDA内存管理和GPU信息查询函数
 */

#pragma once

#include "phase_field_params.h"

/**
 * @brief 打印GPU设备信息
 *
 * 显示所有可用GPU的详细信息
 */
void print_gpu_info();

/**
 * @brief 检查GPU内存是否足够
 *
 * @param required_bytes 需要的字节数
 * @return 内存足够返回true，否则返回false
 */
bool check_gpu_memory(size_t required_bytes);

/**
 * @brief 初始化GPU数组
 *
 * @param arrays GPU数组结构体
 * @param params 模拟参数
 * @return 成功返回true，失败返回false
 */
bool init_gpu_arrays(GPUArrays& arrays, const PhaseFieldParams& params);

/**
 * @brief 释放GPU数组
 *
 * @param arrays GPU数组结构体
 */
void free_gpu_arrays(GPUArrays& arrays);

/**
 * @brief 重置GPU设备
 *
 * 清理GPU状态，释放所有资源
 */
void reset_gpu_device();

/**
 * @brief 设置GPU设备
 *
 * @param device_id 设备ID（从0开始）
 */
void set_gpu_device(int device_id);
