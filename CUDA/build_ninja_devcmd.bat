@echo off
REM ========================================
REM Phase Field CUDA - Build (Ninja + VS DevCmd)
REM ========================================

REM Set UTF-8 code page for Chinese character support
chcp 65001 >nul 2>&1

echo.
echo ========================================
echo Phase Field CUDA - Build (Ninja)
echo ========================================
echo.

REM Set VS path
set "VSPATH=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"

REM Set CUDA path
set "CUDA_PATH=D:\NVIDIA GPU Computing Toolkit\CUDA\v13.1"
set "PATH=%CUDA_PATH%\bin;%PATH%"

echo [0/5] Environment Setup
echo   VS: %VSPATH%
echo   CUDA: %CUDA_PATH%
echo.

REM Call VS Developer Command Prompt
if not exist "%VSPATH%" (
    echo ERROR: VS not found at %VSPATH%
    pause
    exit /b 1
)

call "%VSPATH%"

echo [1/5] Cleaning build directory
if exist build rmdir /s /q build
mkdir build
cd build
echo   [OK] Build directory ready
echo.

echo [2/5] Running CMake with Ninja generator
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86

if errorlevel 1 goto error
echo   [OK] CMake configuration done
echo.

echo [3/5] Compiling with Ninja
ninja

if errorlevel 1 goto error
echo   [OK] Build completed
echo.

echo [4/5] Checking executable
if not exist simulation_cuda.exe (
    echo ERROR: Executable not found
    goto error
)
echo   [OK] Found: simulation_cuda.exe
echo.

echo [5/5] Running simulation
REM Set UTF-8 code page before running the program
chcp 65001 >nul 2>&1
simulation_cuda.exe

echo.
echo ========================================
echo Build and Run completed!
echo ========================================
echo.
pause
exit /b 0

:error
echo.
echo ========================================
echo ERROR: Build failed!
echo ========================================
echo.
pause
exit /b 1
