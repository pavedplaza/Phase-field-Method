/**
 * @file simulation_driver.cpp
 * @brief 相场模拟主程序
 * @author pavedplaza
 * @date 2025-02-15
 *
 * 这是相场法CUDA实现的主入口程序
 * 负责初始化、运行主循环和输出结果
 */

#include "phase_field_params.h"
#include "cuda_headers.h"
#include "cuda_utils.h"
#include "timer.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>  // For std::swap
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _WIN32
#include <windows.h>  // For SetConsoleOutputCP
#include <direct.h>   // For _mkdir
#endif

// CUDA kernel函数声明
extern void launch_init_kernels(GPUArrays& gpu, PhaseFieldParams params);
extern void launch_evolution_kernels(GPUArrays& gpu, PhaseFieldParams params);

// 注释掉Host端初始化函数，现在使用GPU端初始化
/*
void initialize_fields(float* h_phi, float* h_U, const PhaseFieldParams& params) {
    printf("初始化场变量...\n");

    float center_x = params.Nx / 2.0f;
    float center_y = params.Ny / 2.0f;
    float nucleus_radius = 3.2f;  // 晶核半径

    // 使用OpenMP并行初始化
    #pragma omp parallel for collapse(2) if(params.Nx * params.Ny > 10000)
    for (int j = 0; j < params.Ny; j++) {
        for (int i = 0; i < params.Nx; i++) {
            int idx = j * params.Nx + i;

            // 计算到中心的距离
            float x = (i - center_x) * params.dx;
            float y = (j - center_y) * params.dy;
            float r = sqrtf(x * x + y * y);

            // 初始条件：中心晶核
            if (r <= nucleus_radius) {
                h_phi[idx] = 1.0f;  // 固相
            } else {
                // 液相
                h_phi[idx] = -1.0f;
            }

            // 初始浓度场（均匀）
            h_U[idx] = 0.0f;
        }
    }

    printf("  初始条件: 中心晶核 (半径=%.1f)\n", nucleus_radius);
}
*/

/**
 * @brief 主模拟函数
 *
 * 执行相场模拟的主循环
 *
 * @param params 模拟参数
 * @return 成功返回0，失败返回-1
 */
int run_simulation(const PhaseFieldParams& params) {
    printf("\n========================================\n");
    printf("开始相场模拟\n");
    printf("========================================\n");

    // 计算总时间步数
    int num_steps = static_cast<int>(params.total_time / params.dtime);
    printf("总时间步数: %d\n", num_steps);

    // ========== 1. 分配GPU内存并初始化 ==========
    printf("\n[步骤 1/4] 分配GPU内存...\n");
    GPUArrays arrays;
    if (!init_gpu_arrays(arrays, params)) {
        fprintf(stderr, "错误：GPU内存分配失败\n");
        return EXIT_FAILURE;
    }

    // ========== 2. GPU端初始化场变量 ==========
    printf("\n[步骤 2/4] 初始化场变量（GPU）...\n");
    launch_init_kernels(arrays, params);
    printf("  GPU初始化完成：液相φ=-1，中心晶核φ=1\n");

    // ========== 3. 主时间循环 ==========
    printf("\n[步骤 3/4] 开始时间演化...\n");
    printf("────────────────────────────────────────\n");

    Timer total_timer;
    total_timer.start();

    int output_counter = 0;
    int frame_counter = 0;  // 视频帧计数器

    // 文件名缓冲区
    char phi_filename[256];
    char U_filename[256];

    // 创建输出目录
    printf("\n创建输出目录...\n");
    #ifdef _WIN32
    _mkdir("frames");
    #else
    mkdir("frames", 0777);
    #endif
    printf("中间结果将保存在 frames/ 目录\n");

    // 计算保存间隔（每0.5τ0保存一次）
    float save_interval = 0.5f;  // 0.5τ0
    int save_interval_steps = static_cast<int>(save_interval / params.dtime);
    printf("保存间隔: 每%.1fτ0 (约每%d步)\n", save_interval, save_interval_steps);

    for (int step = 0; step < num_steps; step++) {
        // 调用CUDA kernels进行演化
        // 注意：launch_evolution_kernels内部已经完成了指针交换
        launch_evolution_kernels(arrays, params);

        float current_time = step * params.dtime;

        // === 保存中间结果（每0.5τ0）===
        if (step % save_interval_steps == 0 || step == num_steps - 1) {
            // 分配临时Host数组
            size_t total_size = params.Nx * params.Ny * sizeof(float);
            float* h_temp = new float[params.Nx * params.Ny];

            // 保存相场
            CUDA_CHECK(cudaMemcpy(h_temp, arrays.d_phi, total_size, cudaMemcpyDeviceToHost));
            char phi_filename[256];
            sprintf(phi_filename, "frames/phi_frame_%04d.bin", frame_counter);
            FILE* fp_phi = fopen(phi_filename, "wb");
            if (fp_phi) {
                fwrite(h_temp, sizeof(float), params.Nx * params.Ny, fp_phi);
                fclose(fp_phi);
            }

            // 保存浓度场
            CUDA_CHECK(cudaMemcpy(h_temp, arrays.d_U, total_size, cudaMemcpyDeviceToHost));
            char U_filename[256];
            sprintf(U_filename, "frames/U_frame_%04d.bin", frame_counter);
            FILE* fp_U = fopen(U_filename, "wb");
            if (fp_U) {
                fwrite(h_temp, sizeof(float), params.Nx * params.Ny, fp_U);
                fclose(fp_U);
            }

            delete[] h_temp;

            printf("  已保存帧 %4d (step=%d, t=%.3fτ₀)\n",
                   frame_counter, step, current_time);
            frame_counter++;
        }

        // 打印进度
        if (params.enable_output && step % params.output_interval == 0) {
            float progress = 100.0f * step / num_steps;
            printf("  步数 %6d/%6d (%5.1f%%)  时间=%.3f/%.1f\n",
                   step, num_steps, progress, current_time, params.total_time);
            output_counter++;
        }
    }

    total_timer.stop();
    printf("────────────────────────────────────────\n");
    printf("时间演化完成！\n");
    printf("  总耗时: ");
    total_timer.print("总耗时");

    // ========== 4. 结果传回Host并保存 ==========
    printf("\n[步骤 4/4] 传输结果回Host并保存...\n");

    // 分配临时Host数组用于接收结果
    size_t total_size = params.Nx * params.Ny * sizeof(float);
    float* h_phi = new float[params.Nx * params.Ny];
    float* h_U = new float[params.Nx * params.Ny];

    // 从GPU拷贝结果到Host
    CUDA_CHECK(cudaMemcpy(h_phi, arrays.d_phi, total_size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_U, arrays.d_U, total_size, cudaMemcpyDeviceToHost));
    printf("  结果传输完成\n");

    // 保存相场
    FILE* fp_phi = fopen("phi_final.bin", "wb");
    if (fp_phi != nullptr) {
        fwrite(h_phi, sizeof(float), params.Nx * params.Ny, fp_phi);
        fclose(fp_phi);
        printf("  相场已保存: phi_final.bin (%dx%d)\n", params.Nx, params.Ny);
    } else {
        fprintf(stderr, "警告：无法保存相场文件\n");
    }

    // 保存浓度场
    FILE* fp_U = fopen("U_final.bin", "wb");
    if (fp_U != nullptr) {
        fwrite(h_U, sizeof(float), params.Nx * params.Ny, fp_U);
        fclose(fp_U);
        printf("  浓度场已保存: U_final.bin (%dx%d)\n", params.Nx, params.Ny);
    } else {
        fprintf(stderr, "警告：无法保存浓度场文件\n");
    }

    // ========== 5. 清理 ==========
    printf("\n清理资源...\n");

    // 释放Host临时数组
    delete[] h_phi;
    delete[] h_U;

    // 释放GPU内存
    free_gpu_arrays(arrays);

    printf("========================================\n");
    printf("模拟完成！\n");
    printf("========================================\n");

    return 0;
}

/**
 * @brief 主函数
 *
 * 程序入口点
 *
 * @param argc 参数个数
 * @param argv 参数列表
 * @return 成功返回0，失败返回非0
 */
int main(int argc, char* argv[]) {
    // Set console to UTF-8 mode for proper Chinese character display
    #ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
    #endif

    printf("╔════════════════════════════════════════╗\n");
    printf("║   相场法CUDA实现 v1.0                  ║\n");
    printf("║   Phase Field Dendrite Growth         ║\n");
    printf("╚════════════════════════════════════════╝\n\n");

    // 打印GPU信息
    print_gpu_info();

    // 使用默认参数或从配置文件读取
    PhaseFieldParams params;

    // 根据CFL条件计算时间步长
    // dtime = dx^2 / (4 * D) * safety_factor
    params.dtime = params.dx * params.dx / (4.0f * params.D_coefficient) * params.safety_factor;

    // TODO: 添加命令行参数解析
    // TODO: 添加从配置文件读取参数

    printf("模拟参数:\n");
    printf("  网格尺寸:           %dx%d\n", params.Nx, params.Ny);
    printf("  空间步长:           dx=%.3f, dy=%.3f\n", params.dx, params.dy);
    printf("  时间步长:           %.6f (CFL条件: dx²/(4D) × %.2f)\n",
           params.dtime, params.safety_factor);
    printf("  总模拟时间:         %.1f\n", params.total_time);
    printf("  分配系数:           %.3f\n", params.k_partition);
    printf("  耦合系数:           %.2f\n", params.lambda_coup);
    printf("  各向异性强度:       %.3f\n", params.epsilon_aniso);
    printf("  过冷度:             %.2f\n", params.theta_field);
    printf("  周期性边界:         %s\n", params.periodic_boundary ? "是" : "否");
    printf("\n");

    // 检查GPU内存是否足够
    size_t required_memory = 4 * params.Nx * params.Ny * sizeof(float);  // 4个数组
    if (!check_gpu_memory(required_memory)) {
        fprintf(stderr, "\n错误：GPU内存不足，请减小网格尺寸\n");
        return EXIT_FAILURE;
    }

    // 运行模拟
    int result = run_simulation(params);

    if (result == 0) {
        printf("\n✓ 模拟成功完成！\n");
        printf("\n下一步:\n");
        printf("  1. 使用MATLAB读取并可视化结果:\n");
        printf("     load('phi_final.bin');\n");
        printf("     phi = reshape(phi', %d, %d);\n", params.Ny, params.Nx);
        printf("     imagesc(phi'); axis equal tight; colorbar;\n");
        printf("\n  2. 或者使用Python/Matplotlib绘制\n");
        printf("\n");
        return EXIT_SUCCESS;
    } else {
        fprintf(stderr, "\n✗ 模拟失败！\n");
        return EXIT_FAILURE;
    }
}