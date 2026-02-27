/**
 * @file timer.h
 * @brief 高精度计时器类
 * @author pavedplaza
 * @date 2025-02-15
 *
 * 提供高精度的CPU计时功能
 */

#pragma once

#include <chrono>

/**
 * @brief 高精度计时器类
 *
 * 使用C++11的chrono库提供微秒级精度的计时功能
 *
 * 用法示例：
 * @code
 * Timer timer;
 * timer.start();
 * // ... 执行一些操作 ...
 * timer.stop();
 * printf("耗时: %.3f 秒\n", timer.elapsed());
 * @endcode
 */
class Timer {
public:
    /**
     * @brief 构造函数
     */
    Timer();

    /**
     * @brief 析构函数
     */
    ~Timer();

    /**
     * @brief 启动计时器
     */
    void start();

    /**
     * @brief 停止计时器
     */
    void stop();

    /**
     * @brief 重置计时器
     */
    void reset();

    /**
     * @brief 获取经过的时间（秒）
     *
     * @return 经过的时间（浮点秒数）
     */
    double elapsed() const;

    /**
     * @brief 获取经过的时间（毫秒）
     *
     * @return 经过的时间（毫秒）
     */
    double elapsed_milliseconds() const;

    /**
     * @brief 打印经过的时间
     *
     * @param prefix 前缀字符串（用于标识）
     */
    void print(const char* prefix = "耗时") const;

    /**
     * @brief 检查计时器是否正在运行
     *
     * @return 如果正在运行返回true，否则返回false
     */
    bool is_running() const;

private:
    std::chrono::high_resolution_clock::time_point start_time_;  ///< 开始时间点
    std::chrono::high_resolution_clock::time_point end_time_;    ///< 结束时间点
    bool is_running_;                                            ///< 运行状态标志
};
