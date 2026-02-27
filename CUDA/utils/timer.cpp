/**
 * @file timer.cpp
 * @brief 高精度计时器类实现
 * @author pavedplaza
 * @date 2025-02-15
 *
 * 提供高精度的CPU计时功能
 */

#include "timer.h"
#include <cstdio>
#include <cstring>
#include <cmath>

/**
 * @brief 构造函数
 *
 * 初始化计时器状态
 */
Timer::Timer()
    : start_time_()
    , end_time_()
    , is_running_(false)
{
}

/**
 * @brief 析构函数
 */
Timer::~Timer() {
    // 无需特别清理
}

/**
 * @brief 启动计时器
 *
 * 记录开始时间点
 */
void Timer::start() {
    if (is_running_) {
        fprintf(stderr, "警告：计时器已经在运行中\n");
        return;
    }

    start_time_ = std::chrono::high_resolution_clock::now();
    is_running_ = true;
}

/**
 * @brief 停止计时器
 *
 * 记录结束时间点
 */
void Timer::stop() {
    if (!is_running_) {
        fprintf(stderr, "警告：计时器未启动，无法停止\n");
        return;
    }

    end_time_ = std::chrono::high_resolution_clock::now();
    is_running_ = false;
}

/**
 * @brief 重置计时器
 *
 * 清除计时状态
 */
void Timer::reset() {
    start_time_ = std::chrono::high_resolution_clock::time_point{};
    end_time_ = std::chrono::high_resolution_clock::time_point{};
    is_running_ = false;
}

/**
 * @brief 获取经过的时间（秒）
 *
 * @return 经过的时间（浮点秒数）
 */
double Timer::elapsed() const {
    if (is_running_) {
        auto current_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            current_time - start_time_);
        return duration.count() / 1000000.0;
    } else {
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end_time_ - start_time_);
        return duration.count() / 1000000.0;
    }
}

/**
 * @brief 获取经过的时间（毫秒）
 *
 * @return 经过的时间（毫秒）
 */
double Timer::elapsed_milliseconds() const {
    return elapsed() * 1000.0;
}

/**
 * @brief 打印经过的时间
 *
 * 以友好的格式打印时间信息
 *
 * @param prefix 前缀字符串（用于标识）
 */
void Timer::print(const char* prefix) const {
    double time_sec = elapsed();

    if (time_sec < 1.0) {
        printf("%s: %.2f 毫秒\n", prefix, elapsed_milliseconds());
    } else if (time_sec < 60.0) {
        printf("%s: %.3f 秒\n", prefix, time_sec);
    } else if (time_sec < 3600.0) {
        int minutes = static_cast<int>(time_sec / 60.0);
        double seconds = std::fmod(time_sec, 60.0);
        printf("%s: %d 分 %.2f 秒\n", prefix, minutes, seconds);
    } else {
        int hours = static_cast<int>(time_sec / 3600.0);
        int minutes = static_cast<int>(std::fmod(time_sec, 3600.0) / 60.0);
        double seconds = std::fmod(time_sec, 60.0);
        printf("%s: %d 小时 %d 分 %.2f 秒\n", prefix, hours, minutes, seconds);
    }
}

/**
 * @brief 检查计时器是否正在运行
 *
 * @return 如果正在运行返回true，否则返回false
 */
bool Timer::is_running() const {
    return is_running_;
}
