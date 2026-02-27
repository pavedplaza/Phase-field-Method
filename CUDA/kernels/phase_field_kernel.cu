/**
 * @file phase_field_kernel.cu
 * @brief 相场法CUDA Kernels实现
 * @author PhaseFieldCUDA Team
 * @date 2025-02-26
 *
 * Phase 1: 基础CUDA kernel框架
 * - 初始化kernels（完整实现）
 * - 演化kernels（占位实现，完整物理在Phase 2-3实现）
 */

#include "cuda_headers.h"
#include "phase_field_params.h"

// ============================================================
// 初始化Kernels（Phase 1完整实现）
// ============================================================

/**
 * @brief Kernel: 初始化相场数组为-1
 *
 * @param phi 相场数组
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
__global__ void init_phase_field_kernel(float* phi, int Nx, int Ny) {
    // 计算全局线程索引
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // 边界检查
    if (i < Nx && j < Ny) {
        int idx = j * Nx + i;
        phi[idx] = -1.0f;
    }
}

/**
 * @brief Kernel: 初始化浓度场数组为0
 *
 * @param U 浓度场数组
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
__global__ void init_concentration_kernel(float* U, int Nx, int Ny) {
    // 计算全局线程索引
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // 边界检查
    if (i < Nx && j < Ny) {
        int idx = j * Nx + i;
        U[idx] = 0.0f;
    }
}

/**
 * @brief Kernel: 设置中心晶核为1（固相）
 *
 * 在网格中心创建一个圆形晶核区域
 *
 * @param phi 相场数组
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 * @param dx x方向空间步长
 * @param dy y方向空间步长
 * @param nucleus_radius 晶核半径
 */
__global__ void set_nucleus_kernel(float* phi, int Nx, int Ny, float dx, float dy, float nucleus_radius) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < Nx && j < Ny) {
        // 计算到网格中心的距离
        float center_x = Nx * dx / 2.0f;
        float center_y = Ny * dy / 2.0f;
        float x = i * dx;
        float y = j * dy;

        float dist = sqrtf((x - center_x) * (x - center_x) +
                          (y - center_y) * (y - center_y));

        // 如果在晶核半径内，设置为固相
        if (dist <= nucleus_radius) {
            int idx = j * Nx + i;
            phi[idx] = 1.0f;  // 固相
        }
    }
}

// ============================================================
// 设备辅助函数（Phase 2实现）
// ============================================================

/**
 * @brief 设备函数: 周期性边界索引计算
 *
 * 处理负数索引，确保索引在[0, N-1]范围内
 *
 * @param i 原始索引
 * @param N 数组长度
 * @return 处理后的索引
 */
__device__ int periodic_index(int i, int N) {
    return (i + N) % N;  // 处理负数索引
}

/**
 * @brief 设备函数: 计算相场梯度
 *
 * 使用中心差分计算相场在x和y方向的梯度
 *
 * @param phi 相场数组
 * @param i 当前节点的x索引
 * @param j 当前节点的y索引
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 * @param dx x方向空间步长
 * @param dy y方向空间步长
 * @param grad_x 输出: x方向梯度 ∂φ/∂x
 * @param grad_y 输出: y方向梯度 ∂φ/∂y
 */
__device__ void compute_gradient(
    const float* phi, int i, int j, int Nx, int Ny,
    float dx, float dy, float* grad_x, float* grad_y) {

    // 周期性边界索引
    int ip = periodic_index(i + 1, Nx);
    int im = periodic_index(i - 1, Nx);
    int jp = periodic_index(j + 1, Ny);
    int jm = periodic_index(j - 1, Ny);

    // 中心差分
    *grad_x = (phi[j * Nx + ip] - phi[j * Nx + im]) / (2.0f * dx);
    *grad_y = (phi[jp * Nx + i] - phi[jm * Nx + i]) / (2.0f * dy);
}

/**
 * @brief 设备函数: 计算各向异性参数 ε(ψ)
 *
 * ε(ψ) = 1 + ε_aniso · cos[m(ψ - ψ₀)]
 *
 * @param grad_x x方向梯度
 * @param grad_y y方向梯度
 * @param params 模拟参数
 * @return 各向异性参数 ε
 */
__device__ float compute_epsilon(float grad_x, float grad_y, const PhaseFieldParams& params) {
    // ψ = atan2(∂φ/∂y, ∂φ/∂x)
    float angle = atan2f(grad_y, grad_x);

    // ε(ψ) = 1 + ε_aniso · cos[m(ψ - ψ₀)]
    return 1.0f + params.epsilon_aniso * cosf(params.m_aniso * (angle - params.omega0_aniso));
}

/**
 * @brief 设备函数: 计算各向异性参数导数 ε'(ψ)
 *
 * ε'(ψ) = -ε_aniso · m · sin[m(ψ - ψ₀)]
 *
 * @param grad_x x方向梯度
 * @param grad_y y方向梯度
 * @param params 模拟参数
 * @return 各向异性参数导数 ε'
 */
__device__ float compute_epsilon_deriv(float grad_x, float grad_y, const PhaseFieldParams& params) {
    float angle = atan2f(grad_y, grad_x);

    // ε'(ψ) = -ε_aniso · m · sin[m(ψ - ψ₀)]
    return -params.epsilon_aniso * params.m_aniso *
           sinf(params.m_aniso * (angle - params.omega0_aniso));
}

// ============================================================
// Phase 3 设备辅助函数
// ============================================================

/**
 * @brief 设备函数: 计算梯度幅值
 *
 * @param gx x方向梯度
 * @param gy y方向梯度
 * @return 梯度幅值 |∇φ| = sqrt(gx² + gy²)
 */
__device__ float compute_gradient_magnitude(float gx, float gy) {
    return sqrtf(gx * gx + gy * gy);
}

// ============================================================
// Phase 2 Kernels（完整物理实现）
// ============================================================

/**
 * @brief Kernel: 计算相场梯度和各向异性参数
 *
 * 为每个节点计算：
 * - 相场梯度 (∂φ/∂x, ∂φ/∂y)
 * - 各向异性参数 ε(ψ) = 1 + ε·cos[m(ψ-ψ₀)]
 * - 各向异性导数 ε'(ψ) = -ε·m·sin[m(ψ-ψ₀)]
 *
 * @param phi 相场数组
 * @param phidx 输出: x方向梯度数组
 * @param phidy 输出: y方向梯度数组
 * @param epsilon 输出: 各向异性参数数组
 * @param epsilon_deriv 输出: 各向异性导数数组
 * @param params 模拟参数
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
__global__ void compute_gradients_and_epsilon_kernel(
    const float* phi,
    float* phidx, float* phidy,
    float* epsilon, float* epsilon_deriv,
    PhaseFieldParams params, int Nx, int Ny) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < Nx && j < Ny) {
        int idx = j * Nx + i;

        // 计算梯度
        float grad_x, grad_y;
        compute_gradient(phi, i, j, Nx, Ny, params.dx, params.dy, &grad_x, &grad_y);

        phidx[idx] = grad_x;
        phidy[idx] = grad_y;

        // 计算各向异性参数（每个节点根据自己梯度方向计算）
        epsilon[idx] = compute_epsilon(grad_x, grad_y, params);
        epsilon_deriv[idx] = compute_epsilon_deriv(grad_x, grad_y, params);
    }
}

// ============================================================
// 演化Kernels（Phase 1占位实现）
// ============================================================

/**
 * @brief Kernel: 相场演化（Phase 2完整物理实现）
 *
 * 实现完整的相场演化方程（5项）：
 * A² · [k · (1 + (1-k) · U)] · ∂φ/∂t =
 *     ∇·[A² ∇φ]                              ← 项1：扩散项
 *     - ∂/∂x[A·A'·∂φ/∂y]                      ← 项2：各向异性第1部分
 *     + ∂/∂y[A·A'·∂φ/∂x]                      ← 项3：各向异性第2部分
 *     + φ(1 - φ²)                             ← 项4：双阱势项
 *     - λ(1 - φ²)²(θ + kU)                    ← 项5：驱动力项
 *
 * @param phi_new 相场（下一时刻）
 * @param phi 相场（当前时刻）
 * @param U 浓度场
 * @param phidx x方向梯度数组 ∂φ/∂x
 * @param phidy y方向梯度数组 ∂φ/∂y
 * @param epsilon 各向异性参数数组 ε(ψ)
 * @param epsilon_deriv 各向异性导数数组 ε'(ψ)
 * @param params 模拟参数
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
__global__ void evolve_phase_field_kernel(
    float* phi_new, const float* phi, const float* U,
    const float* phidx, const float* phidy,
    const float* epsilon, const float* epsilon_deriv,
    PhaseFieldParams params, int Nx, int Ny) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < Nx && j < Ny) {
        int idx = j * Nx + i;
        int ip = periodic_index(i + 1, Nx);
        int im = periodic_index(i - 1, Nx);
        int jp = periodic_index(j + 1, Ny);
        int jm = periodic_index(j - 1, Ny);

        float phi_center = phi[idx];
        float A_center = epsilon[idx];
        float A2_center = A_center * A_center;

        // === 项1: 扩散项 ∇·[A²∇φ] ===
        float A2_east_avg = 0.5f * (A2_center + epsilon[j * Nx + ip] * epsilon[j * Nx + ip]);
        float A2_west_avg = 0.5f * (A2_center + epsilon[j * Nx + im] * epsilon[j * Nx + im]);
        float A2_north_avg = 0.5f * (A2_center + epsilon[jp * Nx + i] * epsilon[jp * Nx + i]);
        float A2_south_avg = 0.5f * (A2_center + epsilon[jm * Nx + i] * epsilon[jm * Nx + i]);

        float flux_x_east = A2_east_avg * (phi[j * Nx + ip] - phi_center) / params.dx;
        float flux_x_west = A2_west_avg * (phi_center - phi[j * Nx + im]) / params.dx;
        float div_x = (flux_x_east - flux_x_west) / params.dx;

        float flux_y_north = A2_north_avg * (phi[jp * Nx + i] - phi_center) / params.dy;
        float flux_y_south = A2_south_avg * (phi_center - phi[jm * Nx + i]) / params.dy;
        float div_y = (flux_y_north - flux_y_south) / params.dy;

        float diffusion_term = div_x + div_y;

        // === 项2+3: 各向异性项 -∂/∂x[A·A'·∂φ/∂y] + ∂/∂y[A·A'·∂φ/∂x] ===
        float grad_y_east = (phidy[j * Nx + ip] + phidy[j * Nx + i]) * 0.5f;
        float grad_y_west = (phidy[j * Nx + i] + phidy[j * Nx + im]) * 0.5f;
        float AA_grad_y_east = epsilon[j * Nx + ip] * epsilon_deriv[j * Nx + ip] * grad_y_east;
        float AA_grad_y_west = epsilon[j * Nx + im] * epsilon_deriv[j * Nx + im] * grad_y_west;
        float dAdphidy_dx = (AA_grad_y_east - AA_grad_y_west) / (2.0f * params.dx);

        float grad_x_north = (phidx[jp * Nx + i] + phidx[j * Nx + i]) * 0.5f;
        float grad_x_south = (phidx[j * Nx + i] + phidx[jm * Nx + i]) * 0.5f;
        float AA_grad_x_north = epsilon[jp * Nx + i] * epsilon_deriv[jp * Nx + i] * grad_x_north;
        float AA_grad_x_south = epsilon[jm * Nx + i] * epsilon_deriv[jm * Nx + i] * grad_x_south;
        float dAdphidx_dy = (AA_grad_x_north - AA_grad_x_south) / (2.0f * params.dy);

        float anisotropy_term = -dAdphidy_dx + dAdphidx_dy;

        // === 项4: 双阱势项 φ(1 - φ²) ===
        float double_well_term = phi_center * (1.0f - phi_center * phi_center);

        // === 项5: 驱动力项 -λ(1-φ²)²(θ + kU) ===
        float U_val = U[idx];
        float driving_force = -params.lambda_coup * (1.0f - phi_center * phi_center) *
                              (1.0f - phi_center * phi_center) *
                              (params.theta_field + params.k_partition * U_val);

        // === 组合所有项 ===
        float RHS = diffusion_term + anisotropy_term + double_well_term + driving_force;

        // === 时间积分（变系数）===
        float time_coeff = A2_center * (params.k_partition *
                           (1.0f + (1.0f - params.k_partition) * U_val));

        phi_new[idx] = phi_center + params.dtime * RHS / time_coeff;

        // φ ∈ [-1, 1]，无需边界约束
    }
}

// ============================================================
// Phase 3 Kernels（浓度场演化）
// ============================================================

/**
 * @brief Kernel: 计算浓度场梯度
 *
 * 为每个节点计算浓度梯度 ∇U = (∂U/∂x, ∂U/∂y)
 * 使用中心差分和周期性边界条件
 *
 * @param U 浓度场（当前时刻）
 * @param Udx 输出: x方向浓度梯度 ∂U/∂x
 * @param Udy 输出: y方向浓度梯度 ∂U/∂y
 * @param params 模拟参数
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
__global__ void compute_concentration_gradients_kernel(
    const float* U,
    float* Udx, float* Udy,
    PhaseFieldParams params, int Nx, int Ny) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < Nx && j < Ny) {
        int idx = j * Nx + i;

        // 计算梯度（复用compute_gradient函数）
        float grad_x, grad_y;
        compute_gradient(U, i, j, Nx, Ny, params.dx, params.dy, &grad_x, &grad_y);

        Udx[idx] = grad_x;
        Udy[idx] = grad_y;
    }
}

/**
 * @brief Kernel: 计算相场变化率 ∂φ/∂t
 *
 * ∂φ/∂t = (φ_new - φ_old) / dt
 *
 * 注意：必须在swap(d_phi, d_phi_new)之前调用此kernel！
 * 因为需要φ_old（swap前的phi）和φ_new（swap前的phi_new）
 *
 * @param phi_old 相场（旧时刻，swap前的d_phi）
 * @param phi_new 相场（新时刻，swap前的d_phi_new）
 * @param phi_change_rate 输出: 相场变化率 ∂φ/∂t
 * @param dtime 时间步长
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
__global__ void compute_phi_change_rate_kernel(
    const float* phi_old, const float* phi_new,
    float* phi_change_rate, float dtime, int Nx, int Ny) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < Nx && j < Ny) {
        int idx = j * Nx + i;

        // ∂φ/∂t = (φ_new - φ_old) / dt
        phi_change_rate[idx] = (phi_new[idx] - phi_old[idx]) / dtime;
    }
}

/**
 * @brief Kernel: 浓度场演化（Phase 3完整物理实现）
 *
 * 实现完整的浓度场演化方程（3项）：
 * ∂U/∂t = (1/S) × [扩散项 + 耦合项散度 + 源项]
 *
 * 其中 S = [(1+k) - (1-k)φ]/2
 *
 * @param U_new 浓度场（下一时刻）
 * @param phi 相场（当前时刻）
 * @param U 浓度场（当前时刻）
 * @param phidx x方向相场梯度 ∂φ/∂x
 * @param phidy y方向相场梯度 ∂φ/∂y
 * @param Udx x方向浓度梯度 ∂U/∂x
 * @param Udy y方向浓度梯度 ∂U/∂y
 * @param phi_change_rate 相场变化率 ∂φ/∂t
 * @param params 模拟参数
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
__global__ void evolve_concentration_kernel(
    float* U_new, const float* phi, const float* U,
    const float* phidx, const float* phidy,
    const float* Udx, const float* Udy,
    const float* phi_change_rate,
    PhaseFieldParams params, int Nx, int Ny) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < Nx && j < Ny) {
        int idx = j * Nx + i;
        int ip = periodic_index(i + 1, Nx);
        int im = periodic_index(i - 1, Nx);
        int jp = periodic_index(j + 1, Ny);
        int jm = periodic_index(j - 1, Ny);

        float phi_center = phi[idx];
        float U_center = U[idx];

        // === 项1: 扩散项 ∇·[D(1-φ)/2 ∇U] ===
        // 使用通量方法，与相场扩散项实现一致
        float D_center = params.D_coefficient * (1.0f - phi_center) / 2.0f;

        // x方向通量
        float D_east = 0.5f * (D_center + params.D_coefficient * (1.0f - phi[j * Nx + ip]) / 2.0f);
        float D_west = 0.5f * (D_center + params.D_coefficient * (1.0f - phi[j * Nx + im]) / 2.0f);

        float flux_x_east = D_east * (U[j * Nx + ip] - U_center) / params.dx;
        float flux_x_west = D_west * (U_center - U[j * Nx + im]) / params.dx;
        float div_x_diffusion = (flux_x_east - flux_x_west) / params.dx;

        // y方向通量
        float D_north = 0.5f * (D_center + params.D_coefficient * (1.0f - phi[jp * Nx + i]) / 2.0f);
        float D_south = 0.5f * (D_center + params.D_coefficient * (1.0f - phi[jm * Nx + i]) / 2.0f);

        float flux_y_north = D_north * (U[jp * Nx + i] - U_center) / params.dy;
        float flux_y_south = D_south * (U_center - U[jm * Nx + i]) / params.dy;
        float div_y_diffusion = (flux_y_north - flux_y_south) / params.dy;

        float diffusion_term = div_x_diffusion + div_y_diffusion;

        // === 项2: 耦合项散度 ∇·[(1+(1-k)U)/(2√2) (∂φ/∂t)(∇φ/|∇φ|)] ===
        // 需要计算4个方向的梯度幅值
        float grad_mag_east = compute_gradient_magnitude(phidx[j * Nx + ip], phidy[j * Nx + ip]);
        float grad_mag_west = compute_gradient_magnitude(phidx[j * Nx + im], phidy[j * Nx + im]);
        float grad_mag_north = compute_gradient_magnitude(phidx[jp * Nx + i], phidy[jp * Nx + i]);
        float grad_mag_south = compute_gradient_magnitude(phidx[jm * Nx + i], phidy[jm * Nx + i]);

        float coupling_term;

        // 检查梯度幅值，任一方向太小则整个耦合项为0
        if (grad_mag_east < 1e-12f || grad_mag_west < 1e-12f ||
            grad_mag_north < 1e-12f || grad_mag_south < 1e-12f) {
            coupling_term = 0.0f;
        } else {
            // x方向耦合通量
            float coupling_coeff_east = (1.0f + (1.0f - params.k_partition) * U[j * Nx + ip]) / (2.0f * sqrtf(2.0f));
            float coupling_coeff_west = (1.0f + (1.0f - params.k_partition) * U[j * Nx + im]) / (2.0f * sqrtf(2.0f));

            float coupling_flux_x_east = coupling_coeff_east * phi_change_rate[idx] * phidx[j * Nx + ip] / grad_mag_east;
            float coupling_flux_x_west = coupling_coeff_west * phi_change_rate[idx] * phidx[j * Nx + im] / grad_mag_west;

            // y方向耦合通量
            float coupling_coeff_north = (1.0f + (1.0f - params.k_partition) * U[jp * Nx + i]) / (2.0f * sqrtf(2.0f));
            float coupling_coeff_south = (1.0f + (1.0f - params.k_partition) * U[jm * Nx + i]) / (2.0f * sqrtf(2.0f));

            float coupling_flux_y_north = coupling_coeff_north * phi_change_rate[idx] * phidy[jp * Nx + i] / grad_mag_north;
            float coupling_flux_y_south = coupling_coeff_south * phi_change_rate[idx] * phidy[jm * Nx + i] / grad_mag_south;

            // 计算散度
            float div_coupling = (coupling_flux_x_east - coupling_flux_x_west) / (2.0f * params.dx) +
                                 (coupling_flux_y_north - coupling_flux_y_south) / (2.0f * params.dy);
            coupling_term = div_coupling;
        }

        // === 项3: 源项 (1+(1-k)U)/2 × ∂φ/∂t ===
        float source_term = (1.0f + (1.0f - params.k_partition) * U_center) / 2.0f * phi_change_rate[idx];

        // === 项4: 对流项 ===
        // 不考虑速度场，对流项始终为0
        float advection_term = 0.0f;

        // === 组合所有项 ===
        float U_rhs = diffusion_term + coupling_term + source_term + advection_term;

        // === 时间尺度系数 ===
        float S_coefficient = ((1.0f + params.k_partition) - (1.0f - params.k_partition) * phi_center) / 2.0f;

        // === 时间积分：显式欧拉法 ===
        U_new[idx] = U_center + params.dtime * U_rhs / S_coefficient;
    }
}

// ============================================================
// Host包装函数
// ============================================================

/**
 * @brief 启动初始化kernels
 *
 * 初始化GPU上的相场和浓度场数组
 *
 * @param gpu GPU数组结构体
 * @param params 模拟参数
 */
void launch_init_kernels(GPUArrays& gpu, PhaseFieldParams params) {
    // 配置kernel执行参数
    // Block大小: 16x16 = 256线程（适合RTX 3050）
    dim3 block(16, 16);

    // Grid大小: 根据网格尺寸计算，向上取整
    dim3 grid((params.Nx + block.x - 1) / block.x,
              (params.Ny + block.y - 1) / block.y);

    // 启动相场初始化kernel（全部设为液相-1）
    init_phase_field_kernel<<<grid, block>>>(gpu.d_phi, params.Nx, params.Ny);
    CUDA_KERNEL_CHECK();

    // 启动浓度场初始化kernel
    init_concentration_kernel<<<grid, block>>>(gpu.d_U, params.Nx, params.Ny);
    CUDA_KERNEL_CHECK();

    // 设置中心晶核为固相1（晶核半径约3.2个网格单位）
    float nucleus_radius = 3.2f * params.dx;
    set_nucleus_kernel<<<grid, block>>>(gpu.d_phi, params.Nx, params.Ny,
                                        params.dx, params.dy, nucleus_radius);
    CUDA_KERNEL_CHECK();

    // 同步以确保初始化完成
    CUDA_SYNC_CHECK();
}

/**
 * @brief 启动演化kernels（Phase 2+3完整实现）
 *
 * 执行相场和浓度场的时间演化
 * Step 1: 计算相场梯度和epsilon数组
 * Step 2: 更新相场
 * Step 3: 计算相场变化率（必须在swap之前！）
 * Step 4: 交换相场指针
 * Step 5: 计算浓度场梯度
 * Step 6: 更新浓度场（Phase 3完整物理）
 * Step 7: 交换浓度场指针
 *
 * @param gpu GPU数组结构体
 * @param params 模拟参数
 */
void launch_evolution_kernels(GPUArrays& gpu, PhaseFieldParams params) {
    // 配置kernel执行参数
    dim3 block(16, 16);
    dim3 grid((params.Nx + block.x - 1) / block.x,
              (params.Ny + block.y - 1) / block.y);

    // ========== Step 1: 相场演化（Phase 2已实现）==========

    // Step 1.1: 计算相场梯度和epsilon
    compute_gradients_and_epsilon_kernel<<<grid, block>>>(
        gpu.d_phi, gpu.d_phidx, gpu.d_phidy,
        gpu.d_epsilon, gpu.d_epsilon_deriv,
        params, gpu.Nx, gpu.Ny);
    CUDA_KERNEL_CHECK();

    // Step 1.2: 更新相场
    evolve_phase_field_kernel<<<grid, block>>>(
        gpu.d_phi_new, gpu.d_phi, gpu.d_U,
        gpu.d_phidx, gpu.d_phidy,
        gpu.d_epsilon, gpu.d_epsilon_deriv,
        params, gpu.Nx, gpu.Ny);
    CUDA_KERNEL_CHECK();

    // Step 1.3: 计算相场变化率（必须在swap之前！）
    compute_phi_change_rate_kernel<<<grid, block>>>(
        gpu.d_phi, gpu.d_phi_new, gpu.d_phi_change_rate, params.dtime,
        gpu.Nx, gpu.Ny);
    CUDA_KERNEL_CHECK();

    // Step 1.4: 交换相场指针
    std::swap(gpu.d_phi, gpu.d_phi_new);

    // ========== Step 2: 浓度场演化（Phase 3新增）==========

    // Step 2.1: 计算浓度场梯度
    compute_concentration_gradients_kernel<<<grid, block>>>(
        gpu.d_U, gpu.d_Udx, gpu.d_Udy,
        params, gpu.Nx, gpu.Ny);
    CUDA_KERNEL_CHECK();

    // Step 2.2: 计算浓度场演化（3项物理 + 时间积分）
    evolve_concentration_kernel<<<grid, block>>>(
        gpu.d_U_new, gpu.d_phi, gpu.d_U,
        gpu.d_phidx, gpu.d_phidy, gpu.d_Udx, gpu.d_Udy,
        gpu.d_phi_change_rate,
        params, gpu.Nx, gpu.Ny);
    CUDA_KERNEL_CHECK();

    // Step 2.3: 交换浓度场指针
    std::swap(gpu.d_U, gpu.d_U_new);
}

// ============================================================
// 辅助函数
// ============================================================

/**
 * @brief 计算各向异性参数数组
 *
 * Phase 1: 空实现
 * Phase 2: 实现各向异性参数预计算
 *
 * @param gpu GPU数组结构体
 * @param params 模拟参数
 */
void compute_anisotropy_arrays(GPUArrays& gpu, PhaseFieldParams params) {
    // Phase 1: 空实现（Phase 2实现）
    // Phase 2 将预计算每个节点的各向异性参数 A(ψ) 和 A'(ψ)
}

/**
 * @brief 内存拷贝辅助函数: Host -> Device
 *
 * @param d_ptr Device端指针
 * @param h_ptr Host端指针
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
void copy_host_to_device(float* d_ptr, const float* h_ptr, int Nx, int Ny) {
    size_t size = Nx * Ny * sizeof(float);
    CUDA_CHECK(cudaMemcpy(d_ptr, h_ptr, size, cudaMemcpyHostToDevice));
}

/**
 * @brief 内存拷贝辅助函数: Device -> Host
 *
 * @param h_ptr Host端指针
 * @param d_ptr Device端指针
 * @param Nx x方向网格数
 * @param Ny y方向网格数
 */
void copy_device_to_host(float* h_ptr, const float* d_ptr, int Nx, int Ny) {
    size_t size = Nx * Ny * sizeof(float);
    CUDA_CHECK(cudaMemcpy(h_ptr, d_ptr, size, cudaMemcpyDeviceToHost));
}
