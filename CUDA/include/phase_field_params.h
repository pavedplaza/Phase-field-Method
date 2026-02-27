/**
 * @file phase_field_params.h
 * @brief 相场法模拟参数结构体定义
 * @author pavedplaza
 * @date 2025-02-15
 *
 * 定义了相场模拟所需的所有参数结构体
 */

#pragma once

#include <cstdlib>  // for size_t (MSVC compatible)

/**
 * @brief 相场法物理和数值参数
 *
 * 包含所有材料参数、数值参数和边界条件设置
 */
struct PhaseFieldParams {
    // ========== 材料参数 ==========

    float k_partition;          ///< 溶质分配系数 (典型值: 0.15)
    float epsilon_aniso;        ///< 各向异性强度 (典型值: 0.02)
    float m_aniso;              ///< 各向异性模数，六重对称=6
    float omega0_aniso;         ///< 优先生长方向角度（弧度）
    float lambda_coup;          ///< 相场-浓度场耦合强度 (典型值: 10.0)
    float D_coefficient;        ///< 无量纲扩散系数 (典型值: 6.267)
    float theta_field;          ///< 无量纲过冷度 (典型值: -0.2)

    // ========== 数值参数 ==========

    int Nx;                     ///< x方向网格点数
    int Ny;                     ///< y方向网格点数
    float dx;                   ///< x方向空间步长
    float dy;                   ///< y方向空间步长
    float dtime;                ///< 时间步长（根据CFL条件计算）
    float total_time;           ///< 总模拟时间
    float safety_factor;        ///< 时间步长安全系数（默认0.25）

    // ========== 边界条件 ==========

    bool periodic_boundary;     ///< 是否使用周期性边界条件

    // ========== 控制选项 ==========

    bool enable_output;         ///< 是否输出中间结果
    int output_interval;        ///< 输出间隔（时间步数）

    /**
     * @brief 默认构造函数 - 初始化为典型参数值
     */
    PhaseFieldParams()
        : k_partition(0.15f)
        , epsilon_aniso(0.02f)
        , m_aniso(6.0f)
        , omega0_aniso(0.0f)
        , lambda_coup(10.0f)
        , D_coefficient(6.267f)
        , theta_field(-0.2f)
        , Nx(300)
        , Ny(300)
        , dx(0.8f)
        , dy(0.8f)
        , dtime(0.0f)              // 将根据CFL条件计算
        , total_time(20.0f)
        , safety_factor(0.25f)      // 默认安全系数
        , periodic_boundary(true)
        , enable_output(true)
        , output_interval(100)
    {}
};

/**
 * @brief GPU数组管理结构体
 *
 * 管理所有GPU上的数组和内存分配
 */
struct GPUArrays {
    float *d_phi;               ///< 相场（当前时刻）
    float *d_phi_new;           ///< 相场（下一时刻）
    float *d_U;                 ///< 浓度场（当前时刻）
    float *d_U_new;             ///< 浓度场（下一时刻）
    // Phase 2: 新增梯度数组
    float *d_phidx;             ///< x方向梯度 ∂φ/∂x
    float *d_phidy;             ///< y方向梯度 ∂φ/∂y
    float *d_epsilon;           ///< 各向异性参数数组 ε(ψ) = 1 + ε·cos[m(ψ-ψ₀)]
    float *d_epsilon_deriv;     ///< ε的导数 ε'(ψ) = -ε·m·sin[m(ψ-ψ₀)]
    // Phase 3: 浓度场演化所需数组
    float *d_Udx;               ///< x方向浓度梯度 ∂U/∂x
    float *d_Udy;               ///< y方向浓度梯度 ∂U/∂y
    float *d_phi_change_rate;   ///< 相场变化率 ∂φ/∂t

    size_t pitch;               ///< 内存对齐pitch（用于cudaMallocPitch）
    int Nx;                     ///< x方向网格数
    int Ny;                     ///< y方向网格数

    /**
     * @brief 默认构造函数
     */
    GPUArrays()
        : d_phi(nullptr)
        , d_phi_new(nullptr)
        , d_U(nullptr)
        , d_U_new(nullptr)
        , d_phidx(nullptr)
        , d_phidy(nullptr)
        , d_epsilon(nullptr)
        , d_epsilon_deriv(nullptr)
        , d_Udx(nullptr)
        , d_Udy(nullptr)
        , d_phi_change_rate(nullptr)
        , pitch(0)
        , Nx(0)
        , Ny(0)
    {}
};

/**
 * @brief 模拟状态结构体
 *
 * 跟踪模拟的当前状态
 */
struct SimulationState {
    int current_step;           ///< 当前时间步
    float current_time;         ///< 当前模拟时间
    bool is_initialized;        ///< 是否已初始化

    /**
     * @brief 默认构造函数
     */
    SimulationState()
        : current_step(0)
        , current_time(0.0f)
        , is_initialized(false)
    {}
};

/**
 * @brief 性能统计结构体
 *
 * 记录模拟的性能指标
 */
struct PerformanceStats {
    float total_time;           ///< 总模拟时间（秒）
    float kernel_time;          ///< GPU kernel执行时间（秒）
    float data_transfer_time;   ///< 数据传输时间（秒）
    long long memory_allocated;   ///< 分配的GPU内存（字节）

    /**
     * @brief 默认构造函数
     */
    PerformanceStats()
        : total_time(0.0f)
        , kernel_time(0.0f)
        , data_transfer_time(0.0f)
        , memory_allocated(0)
    {}
};
