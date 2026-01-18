%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
%     严格相场-浓度场耦合显式差分代码修正版         %
%          基于文献参数和正确方程实现              %
%                 左右网格计算2dx                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('PFM_core');

% 版本配置参数
config_version = 'v2';  % 可选: 'v1', 'v2'

% 根据版本配置设置参数
switch config_version
    case 'v1'
        enable_resampling = false;  % 关闭边界粒子重采样
    case 'v2'
        enable_resampling = true;   % 开启边界粒子重采样
    otherwise
        enable_resampling = true;  % 默认开启重采样
end

%== 获取程序开始运行的时间
time0=clock();
format long;

%-- 模拟区域参数设置
Nx = 300;           % x方向网格点数
Ny = 300;           % y方向网格点数
NxNy= Nx*Ny;        % 总网格点数

dx = 0.8;          % x方向空间步长
dy = 0.8;          % y方向空间步长

%--- 时间积分参数设置
total_time = 7.5;      % 总模拟时间
print_interval = 0.5;   % 输出时间间隔
max_steps = 50000;      % 防止无限循环的安全上限

%--- 数值精度和材料特性参数设置
% 定义全局数值精度常量
eps_numerical = 1e-12;      % 数值精度常量，用于稳定性检查

%--- 合金参数设置
% epsilon=0.02, k=0.15

k_partition = 0.15;         % 溶质分配系数
epsilon_aniso = 0.02;       % 各向异性强度
m_aniso = 6;                % 各向异性模数
omega0_aniso = 0;         % 优先生长方向角度

% 耦合常数计算
lambda_coup = 10.0;          % 耦合强度参数
a1 = 0.8839;
W_over_d0 = lambda_coup / a1;

% 扩散系数 - 无量纲化
a2 = 0.6267;
D_coefficient = a2 * lambda_coup;        % 无量纲扩散系数

% 温度场参数 - 基于文献的无量纲过冷度
theta_field = -0.05;         % 无量纲过冷度 (负值促进生长)

% 界面速度 - 对于静态枝晶生长设为0
v_interface = 0;            % 界面速度参数

fprintf('=== 基于文献的严格相场-浓度场耦合模拟 ===\n');
fprintf('网格参数: %dx%d, 区域: %.2fx%.2f\n', Nx, Ny, Nx*dx, Ny*dy);
fprintf('文献参数 (Al-Cu合金):\n');
fprintf('  - 溶质分配系数 k: %.3f\n', k_partition);
fprintf('  - 各向异性强度 ε: %.3f\n', epsilon_aniso);
fprintf('  - 各向异性模数 m: %d (六重对称)\n', m_aniso);
fprintf('  - 优先生长方向 ψ₀: %.1f° (沿x轴)\n', omega0_aniso);
fprintf('  - 耦合强度 λ: %.1f\n', lambda_coup);
fprintf('  - 扩散系数 D: %.3f\n', D_coefficient);
fprintf('  - 无量纲过冷度 θ: %.1f\n', theta_field);

%--- 基于CFL条件计算时间步长
dtime = dx^2 / (4 * D_coefficient);
safety_factor = 0.25;  % 更保守的安全系数
dtime = safety_factor * dtime;

%--- 估计总步数（用于数组预分配）
estimated_max_steps = ceil(total_time / dtime) + 100;

fprintf('时间步长: %.6e\n', dtime);
fprintf('估计总步数: %d\n', estimated_max_steps);

%--- 可视化参数设置
visualization_mode = 1;    % 可视化模式选择
save_video = false;          % 是否保存视频

if save_video
    video_filename = sprintf('corrected_phase_field_%s.mp4', datestr(now, 'yyyymmdd_HHMMSS'));
end

%--- 边界条件设置
boundary_condition_type = 'periodic';    % 'periodic' 或 'zero_gradient'
fprintf('边界条件类型: %s\n', boundary_condition_type);

%--- 数据探测线分析参数设置
enable_probe_analysis = true;     % 启用探测线分析 (开启)
probe_line_position = 0.5;        % 探测线y位置比例(0.5=Ny/2) (开启)
probe_display_mode = 'overlay';  % 显示模式：'overlay'(叠加) 或 'separate'(分离)
save_probe_results = true;        % 保存探测线数据 (开启)
probe_image_format = 'png';       % 探测线图像保存格式 (开启)

%-- 初始化相场和浓度场
fprintf('\n初始化相场和浓度场...\n');

phi = -ones(Ny, Nx);       % 初始化为全液相 (-1)
U = zeros(Ny, Nx);         % 初始均匀浓度场

% 创建初始晶核
seed_radius = 3.2;
cx = round(Nx / 2); 
cy = round(Ny / 2);

% 在中心创建圆形晶核
for i = 1:Nx
    for j = 1:Ny
        if ((i - cx)^2 + (j - cy)^2) < seed_radius^2
            phi(i, j) = 1.0;  % 晶核内部设为固相
        end
    end
end

%--- 定义六方向可视化参考线（仅用于显示，不影响phi场计算）
fprintf('定义六方向生长参考线...\n');

% 定义六个优先生长方向角度（弧度）：30°, 90°, 150°, 210°, 270°, 330°
six_tip_angles = [pi/6, pi/2, 5*pi/6, 7*pi/6, 9*pi/6, 11*pi/6];
six_tip_directions = [cos(six_tip_angles); sin(six_tip_angles)]';  % 6×2方向矩阵

center_x = cx * dx;  % 转换为物理坐标
center_y = cy * dy;

% 创建六方向参考线数据用于可视化
six_reference_lines = cell(6, 1);
for tip_id = 1:6
    angle = six_tip_angles(tip_id);

    % 从中心向外延伸的参考线，长度为计算域的一半
    line_length = min(Nx, Ny) * dx / 2.5;

    % 生成参考线坐标点
    num_points = 100;
    r_values = linspace(0, line_length, num_points);
    line_x = center_x + r_values .* cos(angle);
    line_y = center_y + r_values .* sin(angle);

    % 存储参考线坐标
    six_reference_lines{tip_id} = [line_x', line_y'];
end

phi_init = phi;
U_init = U;

fprintf('初始条件:\n');
fprintf('  - 固相像素数: %d / %d (%.1f%%)\n', ...
        sum(phi(:) > 0), Nx*Ny, 100*sum(phi(:) > 0)/(Nx*Ny));
fprintf('  - 六方向参考线：角度 [0°, 60°, 120°, 180°, 240°, 300°]\n');

%--- 初始化可视化
if visualization_mode >= 1
    initialize_corrected_visualization(Nx, Ny);
    fprintf('可视化已启用\n');

    if save_video
        video_handle = video_manager('init', video_filename, gcf);
    end
end

%--- 初始化界面追踪可视化 (独立窗口)
if visualization_mode >= 1
    initialize_interface_visualization(Nx, Ny, dx, dy);
    fprintf('界面追踪可视化已启用 (figure(2))\n');
end

%--- 预分配历史数据存储（基于时间控制）
max_saves = ceil(total_time / print_interval) + 10;

phi_history = cell(max_saves, 1);
U_history = cell(max_saves, 1);
time_history = zeros(max_saves, 1);

%--- 界面追踪变量初始化
interface_history = zeros(max_saves, 1);      % 存储界面长度演化
contour_history = cell(max_saves, 1);         % 存储界面坐标数据

%--- 六重对称枝晶尖端动力学参数存储初始化
% 六尖端数据存储（第1维：时间，第2维：6个尖端，第3维：参数）
tip_position_history_6 = zeros(max_saves, 6, 2);  % [时间, 尖端ID, 坐标(x,y)]
tip_distance_history_6 = zeros(max_saves, 6);      % [时间, 尖端ID] 距离
tip_velocity_history_6 = zeros(max_saves, 6);      % [时间, 尖端ID] 速度
tip_curvature_history_6 = zeros(max_saves, 6);     % [时间, 尖端ID] 曲率半径
tip_history_buffer_6 = [];                        % 六尖端历史缓冲区结构体数组


% 存储初始状态
phi_history{1} = phi;
U_history{1} = U;
time_history(1) = 0;

%--- 初始界面检测
fprintf('检测初始固液界面...\n');
% 初始界面追踪等值线提取
contour_level = 0;
[Ny, Nx] = size(phi);

% 创建坐标网格
x_coords = ((0:Nx-1) + 0.5) * dx;  % 使用网格中心坐标
y_coords = ((0:Ny-1) + 0.5) * dy;  % 使用网格中心坐标

% 提取指定等值线，使用更安全的方法
[C, h] = contour(x_coords, y_coords, phi, [contour_level contour_level]);

% 如果有等值线数据
if ~isempty(C)
    % 从轮廓矩阵C中提取坐标
    contour_data = [];
    i = 1;
    while i < size(C, 2)
        level = C(1, i);
        num_points = C(2, i);

        if num_points > 0
            x_coords_contour = C(1, i+1:i+num_points);
            y_coords_contour = C(2, i+1:i+num_points);

            if isempty(contour_data)
                contour_data = [x_coords_contour(:), y_coords_contour(:)];
            else
                contour_data = [contour_data; x_coords_contour(:), y_coords_contour(:)];
            end
        end

        i = i + num_points + 1;
    end

    % 筛选有效的等值线点 (避免边界异常)
    valid_points = true(size(contour_data, 1), 1);

    for i = 2:size(contour_data, 1)-1
        % 计算相邻点的距离
        dist_x = contour_data(i, 1) - contour_data(i-1, 1);
        dist_y = contour_data(i, 2) - contour_data(i-1, 2);
        dist = sqrt(dist_x^2 + dist_y^2);

        % 如果距离超过合理范围，标记为无效
        max_reasonable_dist = 3 * max(dx, dy);  % 最大合理距离
        if dist > max_reasonable_dist
            % 检查是否接近边界 (过滤边界异常点)
            if contour_data(i, 1) < 2*dx || contour_data(i, 1) > (Nx-2)*dx || ...
               contour_data(i, 2) < 2*dy || contour_data(i, 2) > (Ny-2)*dy
                valid_points(i) = false;
            end
        end
    end

    % 提取有效点
    if sum(valid_points) > 2  % 至少保留3个点
        contour_data = contour_data(valid_points, :);
    end
else
    contour_data = [];
end

contour_lines = C;  % 返回轮廓矩阵
interface_length = size(contour_data, 1);
interface_history(1) = interface_length;
contour_history{1} = contour_data;
fprintf('初始界面长度: %.2f 网格单位\n', interface_length);

% 初始界面可视化更新
if visualization_mode >= 1
    update_interface_visualization(phi, contour_data, interface_length, 0, time_history, interface_history(1), Nx, Ny, dx, dy, [], [], [], []);
end

%--- 开始演化计算
fprintf('\n=== 基于时间的严格耦合模拟 ===\n');
fprintf('目标模拟时间: %.1f 秒\n', total_time);
fprintf('时间步长: %.6e 秒\n', dtime);
fprintf('估计总步数: %d\n', estimated_max_steps);
fprintf('输出间隔: %.1f 秒\n', print_interval);

% 初始化时间和循环变量
current_time = 0;
istep = 0;
next_print_time = print_interval;
save_count = 1;  % 初始状态已保存在第1个位置

while current_time < total_time && istep < max_steps
    istep = istep + 1;

    % 预分配数组
    lap_phi = zeros(Ny, Nx);
    lap_U = zeros(Ny, Nx);
    phidx = zeros(Ny, Nx);
    phidy = zeros(Ny, Nx);
    Udx = zeros(Ny, Nx);
    Udy = zeros(Ny, Nx);
    epsilon_arr = zeros(Ny, Nx);
    epsilon_deriv_arr = zeros(Ny, Nx);
    interface_mask = false(Ny, Nx);  % 界面标记

    % 第一步：计算梯度和拉普拉斯算子，并检测界面区域
    for i = 1:Nx
        for j = 1:Ny
            % 应用周期性边界条件
            [ip, im, jp, jm] = apply_boundary_conditions(i, j, Nx, Ny , boundary_condition_type);

            % 计算相场拉普拉斯算子
            lap_phi(i,j) = (phi(ip,j) + phi(im,j) + phi(i,jp) + phi(i,jm) - 4*phi(i,j)) / (dx*dx);

            % 计算相场梯度
            phidx_temp = (phi(ip,j) - phi(im,j)) / (2*dx);
            phidy_temp = (phi(i,jp) - phi(i,jm)) / (2*dy);
            phidx(i,j) = phidx_temp;
            phidy(i,j) = phidy_temp;

            % 检测界面：|∇φ| > eps_numerical
            grad_phi_mag = sqrt(phidx_temp^2 + phidy_temp^2);
            if grad_phi_mag > eps_numerical
                interface_mask(i,j) = true;
            end
        end
    end

    % 统计界面区域信息
    interface_points = sum(interface_mask(:));
    interface_percentage = 100.0 * interface_points / (Nx * Ny);

    % 如果没有界面区域，跳过本次计算
    if ~any(interface_mask(:))
        fprintf('时间: %8.4f/%.1f秒, 步数: %4d, 无界面变化，跳过计算\n', current_time, total_time, istep);
        continue;
    end
    
    % 第一个循环：计算剩余空间导数（只在界面区域）
    for i = 1:Nx
        for j = 1:Ny
            % 只计算界面区域的导数
            if ~interface_mask(i,j)
                continue;
            end

            % 应用周期性边界条件
            [ip, im, jp, jm] = apply_boundary_conditions(i, j, Nx, Ny , boundary_condition_type);

            % 计算U场拉普拉斯算子（相场的已在前面计算）
            lap_U(i,j) = (U(ip,j) + U(im,j) + U(i,jp) + U(i,jm) - 4*U(i,j)) / (dx*dx);

            % 计算U场梯度
            Udx(i,j) = (U(ip,j) - U(im,j)) / (2*dx);
            Udy(i,j) = (U(i,jp) - U(i,jm)) / (2*dy);

            % 计算各向异性参数（基于已计算的相场梯度）
            theta = atan2(phidy(i,j), phidx(i,j));
            epsilon_arr(i,j) = 1.0 + epsilon_aniso * cos(m_aniso * (theta - omega0_aniso));
            epsilon_deriv_arr(i,j) = -epsilon_aniso * m_aniso * sin(m_aniso * (theta - omega0_aniso));
        end
    end
    
    % 保存当前phi用于计算时间导数
    phi_old = phi;
    
    % 第二个循环：更新场变量（只在界面区域）
    for i = 1:Nx
        for j = 1:Ny
            % 只更新界面区域的场变量
            if ~interface_mask(i,j)
                continue;
            end

            % 应用周期性边界条件
            [ip, im, jp, jm] = apply_boundary_conditions(i, j, Nx, Ny , boundary_condition_type);
            
            phi_old_val = phi_old(i,j);
            U_old_val = U(i,j);
            epsilon_val = epsilon_arr(i,j);
            epsilon_deriv_val = epsilon_deriv_arr(i,j);

            % === 相场方程计算 ===
            % 检查相场梯度幅值
            grad_phi_mag = sqrt(phidx(i,j)^2 + phidy(i,j)^2);

            if grad_phi_mag < eps_numerical
                % 如果相场梯度太小，各向异性项和扩散项都设为0
                anisotropy_terms = 0;
                diffusion_term = 0;
            else
                % 正常计算各向异性项和扩散项

                % 正确实现各向异性项: -∂/∂x[A(ψ)A'(ψ)∂ϕ/∂y] + ∂/∂y[A(ψ)A'(ψ)∂ϕ/∂x]

                % 计算 ∂/∂x[A(ψ)A'(ψ)∂ϕ/∂y]
                term_A_east = epsilon_arr(ip,j) * epsilon_deriv_arr(ip,j) * phidy(ip,j);
                term_A_west = epsilon_arr(im,j) * epsilon_deriv_arr(im,j) * phidy(im,j);
                dAdphidy_dx = (term_A_east - term_A_west) / (2*dx);

                % 计算 ∂/∂y[A(ψ)A'(ψ)∂ϕ/∂x]
                term_B_north = epsilon_arr(i,jp) * epsilon_deriv_arr(i,jp) * phidx(i,jp);
                term_B_south = epsilon_arr(i,jm) * epsilon_deriv_arr(i,jm) * phidx(i,jm);
                dAdphidx_dy = (term_B_north - term_B_south) / (2*dy);

                anisotropy_terms = -dAdphidy_dx + dAdphidx_dy;

                % 正确计算扩散项: ∇·(A(ψ)²∇ϕ)
                % 使用通量方法计算散度，确保各向异性正确实现
                A2_center = epsilon_val^2;

                % x方向通量：使用界面处的平均A²值（算术平均）
                A2_east = 0.5 * (A2_center + epsilon_arr(ip,j)^2);
                A2_west = 0.5 * (A2_center + epsilon_arr(im,j)^2);

                flux_x_east = A2_east * (phi(ip,j) - phi_old_val) / dx;
                flux_x_west = A2_west * (phi_old_val - phi(im,j)) / dx;
                div_x = (flux_x_east - flux_x_west) / dx;

                % y方向通量
                A2_north = 0.5 * (A2_center + epsilon_arr(i,jp)^2);
                A2_south = 0.5 * (A2_center + epsilon_arr(i,jm)^2);

                flux_y_north = A2_north * (phi(i,jp) - phi_old_val) / dy;
                flux_y_south = A2_south * (phi_old_val - phi(i,jm)) / dy;
                div_y = (flux_y_north - flux_y_south) / dy;

                diffusion_term = div_x + div_y;
            end
            
            % 双阱势项和驱动力项
            double_well_term = phi_old_val * (1.0 - phi_old_val^2);
            driving_force = -lambda_coup * (1.0 - phi_old_val^2)^2 * (theta_field + k_partition * U_old_val);
            
            % 相场右端项
            phi_rhs = diffusion_term + anisotropy_terms + double_well_term + driving_force;
            
            % 时间尺度系数 - 根据理论方程 A(ψ)²[k(1+(1-k)U)] ∂ϕ/∂t
            time_coeff_phi = A2_center * (k_partition * (1.0 + (1.0 - k_partition) * U_old_val));
            % time_coeff_phi =  (k_partition * (1.0 + (1.0 - k_partition) * U_old_val));
            
%             % 相场时间积分
%             if abs(time_coeff_phi) > eps_numerical
%                 phi(i,j) = phi_old_val + dtime * phi_rhs / time_coeff_phi;
%             else
%                 phi(i,j) = phi_old_val;
%             end

            phi(i,j) = phi_old_val + dtime * phi_rhs / time_coeff_phi;

            % === 浓度场方程计算 ===
            phi_change_rate = (phi(i,j) - phi_old_val) / dtime;

            % 项1: 扩散项 ∇·[D(1-ϕ)/2 ∇U]
            % 使用通量方法计算变系数扩散项，参照相场方程计算方法
            D_center = D_coefficient * (1.0 - phi_old_val) / 2.0;

            % x方向通量：使用界面处的平均扩散系数值
            D_east = 0.5 * (D_center + D_coefficient * (1.0 - phi(ip,j)) / 2.0);
            D_west = 0.5 * (D_center + D_coefficient * (1.0 - phi(im,j)) / 2.0);

            flux_x_east = D_east * (U(ip,j) - U_old_val) / dx;
            flux_x_west = D_west * (U_old_val - U(im,j)) / dx;
            div_x_diffusion = (flux_x_east - flux_x_west) / dx;

            % y方向通量
            D_north = 0.5 * (D_center + D_coefficient * (1.0 - phi(i,jp)) / 2.0);
            D_south = 0.5 * (D_center + D_coefficient * (1.0 - phi(i,jm)) / 2.0);

            flux_y_north = D_north * (U(i,jp) - U_old_val) / dy;
            flux_y_south = D_south * (U_old_val - U(i,jm)) / dy;
            div_y_diffusion = (flux_y_north - flux_y_south) / dy;

            div_diffusion = div_x_diffusion + div_y_diffusion;
            
           % 项2: 耦合项散度 ∇·[ (1+(1-k)U)/(2√2) (∂ϕ/∂t)(∇ϕ/|∇ϕ|) ]

            % 修正耦合通量散度计算：使用方向特定的耦合系数
            % 计算各个方向的梯度幅值
            grad_mag_east = sqrt(phidx(ip,j)^2 + phidy(ip,j)^2);
            grad_mag_west = sqrt(phidx(im,j)^2 + phidy(im,j)^2);
            grad_mag_north = sqrt(phidx(i,jp)^2 + phidy(i,jp)^2);
            grad_mag_south = sqrt(phidx(i,jm)^2 + phidy(i,jm)^2);

            % 检查所有方向的梯度幅值，如果任何一个小于eps_numerical，则整个耦合项为0
            if grad_mag_east < eps_numerical || grad_mag_west < eps_numerical || ...
               grad_mag_north < eps_numerical || grad_mag_south < eps_numerical
                % 如果任何方向梯度太小，整个耦合项散度设为0
                div_coupling = 0;
            else
                % 所有方向梯度都足够大，正常计算耦合项散度

                % x方向耦合通量 - 使用对应节点的U值计算耦合系数
                coupling_coeff_east = (1.0 + (1.0 - k_partition) * U(ip,j)) / (2.0 * sqrt(2.0));
                coupling_coeff_west = (1.0 + (1.0 - k_partition) * U(im,j)) / (2.0 * sqrt(2.0));

                coupling_flux_x_east = coupling_coeff_east * phi_change_rate * phidx(ip,j) / grad_mag_east;
                coupling_flux_x_west = coupling_coeff_west * phi_change_rate * phidx(im,j) / grad_mag_west;

                % y方向耦合通量 - 使用对应节点的U值计算耦合系数
                coupling_coeff_north = (1.0 + (1.0 - k_partition) * U(i,jp)) / (2.0 * sqrt(2.0));
                coupling_coeff_south = (1.0 + (1.0 - k_partition) * U(i,jm)) / (2.0 * sqrt(2.0));

                coupling_flux_y_north = coupling_coeff_north * phi_change_rate * phidy(i,jp) / grad_mag_north;
                coupling_flux_y_south = coupling_coeff_south * phi_change_rate * phidy(i,jm) / grad_mag_south;

                % 计算散度
                div_coupling = (coupling_flux_x_east - coupling_flux_x_west)/(2*dx) + ...
                               (coupling_flux_y_north - coupling_flux_y_south)/(2*dy);
            end
            
            % 项3: 源项 (1+(1-k)U)/2 ∂ϕ/∂t
            source_term = (1.0 + (1.0 - k_partition) * U_old_val) / 2.0 * phi_change_rate;
            
            % 项4: 对流通量项 (简化实现，因v_interface=0)
            advection_term = 0; % 当v_interface=0时此项为0
            
            % 浓度场右端项
            U_rhs = div_diffusion + div_coupling + source_term + advection_term;
            
            % 时间尺度系数
            S_coefficient = ((1.0 + k_partition) - (1.0 - k_partition) * phi_old_val) / 2.0;
            
            % 浓度场时间积分
%             if abs(S_coefficient) > eps_numerical
%                 U_new = U_old_val + dtime * U_rhs / S_coefficient;
%                 U(i,j) = U_new;
%             else
%                 U(i,j) = U_old_val;
%             end

            U_new = U_old_val + dtime * U_rhs / S_coefficient;
            U(i,j) = U_new;

        end
    end
    
    % 更新模拟时间
    current_time = current_time + dtime;

    % 处理最后的非完整时间步
    if current_time + dtime > total_time
        final_dtime = total_time - current_time;
        if final_dtime > 1e-12  % 如果还有时间需要处理
            % 使用调整的时间步长进行最后一步
            % 注意：这里简化处理，为了精确到达total_time
            % 在实际应用中，可以实施更精确的部分时间步处理
            current_time = total_time;
        end
    end

    % 输出结果（基于时间间隔）
    if (current_time >= next_print_time || current_time >= total_time) && istep >= 20
        fprintf('时间: %.4f/%.1fs (%.1f%%), 步数: %d, 界面点数: %d (%.2f%%)\n', ...
                current_time, total_time, 100*current_time/total_time, istep, interface_points, interface_percentage);

%--- 界面追踪
% 界面追踪等值线提取
contour_level = 0;
[Ny, Nx] = size(phi);

% 创建坐标网格
x_coords = ((0:Nx-1) + 0.5) * dx;  % 使用网格中心坐标
y_coords = ((0:Ny-1) + 0.5) * dy;  % 使用网格中心坐标

% 提取指定等值线，使用更安全的方法
[C, h] = contour(x_coords, y_coords, phi, [contour_level contour_level]);

% 如果有等值线数据
if ~isempty(C)
    % 从轮廓矩阵C中提取坐标
    contour_data = [];
    i = 1;
    while i < size(C, 2)
        level = C(1, i);
        num_points = C(2, i);

        if num_points > 0
            x_coords_contour = C(1, i+1:i+num_points);
            y_coords_contour = C(2, i+1:i+num_points);

            if isempty(contour_data)
                contour_data = [x_coords_contour(:), y_coords_contour(:)];
            else
                contour_data = [contour_data; x_coords_contour(:), y_coords_contour(:)];
            end
        end

        i = i + num_points + 1;
    end

    % 筛选有效的等值线点 (避免边界异常)
    valid_points = true(size(contour_data, 1), 1);

    for i = 2:size(contour_data, 1)-1
        % 计算相邻点的距离
        dist_x = contour_data(i, 1) - contour_data(i-1, 1);
        dist_y = contour_data(i, 2) - contour_data(i-1, 2);
        dist = sqrt(dist_x^2 + dist_y^2);

        % 如果距离超过合理范围，标记为无效
        max_reasonable_dist = 3 * max(dx, dy);  % 最大合理距离
        if dist > max_reasonable_dist
            % 检查是否接近边界 (过滤边界异常点)
            if contour_data(i, 1) < 2*dx || contour_data(i, 1) > (Nx-2)*dx || ...
               contour_data(i, 2) < 2*dy || contour_data(i, 2) > (Ny-2)*dy
                valid_points(i) = false;
            end
        end
    end

    % 提取有效点
    if sum(valid_points) > 2  % 至少保留3个点
        contour_data = contour_data(valid_points, :);
    end
else
    contour_data = [];
end

contour_lines = C;  % 返回轮廓矩阵
        interface_length = size(contour_data, 1);
        fprintf('界面长度: %.2f 网格单位\n', interface_length);

        % 保存历史数据
        save_count = save_count + 1;
        if save_count <= max_saves
            phi_history{save_count} = phi;
            U_history{save_count} = U;
            time_history(save_count) = current_time;
            interface_history(save_count) = interface_length;
            contour_history{save_count} = contour_data;

            % 更新下次输出时间
            next_print_time = next_print_time + print_interval;

            %--- 六重对称枝晶尖端追踪分析
            fprintf('开始六尖端追踪分析...\n');

            % 边界粒子重采样处理
            if enable_resampling && ~isempty(contour_data) && size(contour_data, 1) > 3
                target_spacing = 0.4;
                contour_data = resample_contour_equidistant(contour_data, target_spacing);
            end

            % 使用六尖端检测算法
            six_tips_info = find_six_tips_direction_specific(...
                contour_data, cx, cy, dx, dy, tip_history_buffer_6, six_tip_angles);

            % 计算六尖端动力学参数
            [six_tips_dynamics, tip_history_buffer_6] = calculate_six_tip_dynamics(...
                six_tips_info, contour_data, tip_history_buffer_6, istep, dtime, cx, cy, dx, dy);

            % 存储六尖端数据
            for tip_id = 1:6
                if six_tips_dynamics(tip_id).is_valid
                    tip_position_history_6(save_count, tip_id, :) = six_tips_dynamics(tip_id).position;
                    tip_distance_history_6(save_count, tip_id) = six_tips_dynamics(tip_id).distance;
                    tip_velocity_history_6(save_count, tip_id) = six_tips_dynamics(tip_id).velocity;
                    tip_curvature_history_6(save_count, tip_id) = six_tips_dynamics(tip_id).curvature_radius;
                else
                    tip_position_history_6(save_count, tip_id, :) = [0, 0];
                    tip_distance_history_6(save_count, tip_id) = 0;
                    tip_velocity_history_6(save_count, tip_id) = 0;
                    tip_curvature_history_6(save_count, tip_id) = 0;
                end
            end
        end

        % 可视化更新
        if visualization_mode >= 1
            % 更新figure(1)：相场、浓度场、C/C∞
            update_corrected_visualization(phi, U, istep, [], Nx, Ny, probe_line_position, k_partition);

            % 更新figure(2)：界面追踪分析
            update_interface_visualization(phi, contour_data, interface_length, istep, time_history, interface_history(1:save_count), Nx, Ny, dx, dy, [], [], [], []);

            if save_video
                video_manager('add', '', gcf);
            end
        end
    end
end

%--- 数组清理（裁剪到实际使用大小）
fprintf('\n清理数组到实际使用大小...\n');
if save_count < max_saves
    phi_history = phi_history(1:save_count);
    U_history = U_history(1:save_count);
    time_history = time_history(1:save_count);
    interface_history = interface_history(1:save_count);
    contour_history = contour_history(1:save_count);
    tip_position_history_6 = tip_position_history_6(1:save_count, :, :);
    tip_distance_history_6 = tip_distance_history_6(1:save_count, :);
    tip_velocity_history_6 = tip_velocity_history_6(1:save_count, :);
    tip_curvature_history_6 = tip_curvature_history_6(1:save_count, :);
    fprintf('数组已裁剪: %d -> %d 个保存点\n', max_saves, save_count);
else
    fprintf('数组大小无需调整: %d 个保存点\n', save_count);
end

fprintf('模拟完成! 实际运行时间: %.4f/%.1f 秒\n', current_time, total_time);

%--- 完成视频录制
if save_video && visualization_mode >= 1
    video_manager('close', '', '');
    fprintf('视频已保存: %s\n', video_filename);
end

%--- 数据探测线分析 (已开启)
if enable_probe_analysis
    fprintf('\n开始数据探测线分析...\n');

    % 提取探测线数据
    [y_indices, CCinf_data, phi_line, U_line] = extract_probe_data(...
        phi, U, k_partition, Nx, Ny, probe_line_position);

    % 计算实际探测线位置
    y_probe = round(Ny * probe_line_position);

    % 根据显示模式选择可视化方式
    switch probe_display_mode
        case 'overlay'
            % 叠加模式：在现有图上添加探测线和独立分析窗口
            plot_probe_analysis_overlay(phi, U, CCinf_data, istep, Nx, Ny, probe_line_position, k_partition);
        case 'separate'
            % 分离模式：创建独立的综合分析窗口
            plot_probe_analysis_separate(phi, U, CCinf_data, istep, Nx, Ny, probe_line_position, k_partition);
        otherwise
            fprintf('未知的探测线显示模式: %s\n', probe_display_mode);
    end

    % 保存探测线数据
    if save_probe_results
        save_probe_analysis_data(y_indices, CCinf_data, phi_line, U_line, istep, y_probe, W_over_d0);
    end

    fprintf('探测线分析完成，位置: y=%d (第%.1f%%高度)\n', y_probe, probe_line_position * 100);
    fprintf('C/C∞范围: [%.3f, %.3f]\n', min(CCinf_data), max(CCinf_data));
end


%--- 六重对称枝晶尖端分析和可视化
fprintf('\n进行六尖端动力学分析...\n');
plot_six_tip_dynamics_analysis(time_history, tip_position_history_6, tip_distance_history_6, tip_velocity_history_6, tip_curvature_history_6, print_interval);

%--- 尖端4专项分析
fprintf('\n进行尖端4专项分析...\n');
% 提取尖端4的有效数据
valid_time_steps = time_history > 0;
valid_times = time_history(valid_time_steps);
if length(valid_times) >= 2 && size(tip_velocity_history_6, 2) >= 4 && size(tip_curvature_history_6, 2) >= 4
    tip4_velocities = tip_velocity_history_6(valid_time_steps, 4);
    tip4_curvatures = tip_curvature_history_6(valid_time_steps, 4);

    % 过滤有效数据
    valid_velocity_mask = abs(tip4_velocities) > 1e-6;
    valid_curvature_mask = isfinite(tip4_curvatures) & tip4_curvatures > 0 & tip4_curvatures < 100;

    if sum(valid_velocity_mask) > 0 || sum(valid_curvature_mask) > 0
        generate_tip_detailed_analysis(valid_times, tip4_velocities, tip4_curvatures, 4);
    else
        fprintf('尖端4数据不足，跳过专项分析\n');
    end
else
    fprintf('尖端4数据不足，跳过专项分析\n');
end

%--- 尖端1专项分析
fprintf('\n进行尖端1专项分析...\n');
if length(valid_times) >= 2 && size(tip_velocity_history_6, 2) >= 1 && size(tip_curvature_history_6, 2) >= 1
    tip1_velocities = tip_velocity_history_6(valid_time_steps, 1);
    tip1_curvatures = tip_curvature_history_6(valid_time_steps, 1);

    % 数据验证
    valid_velocity_mask = abs(tip1_velocities) > 1e-6;
    valid_curvature_mask = isfinite(tip1_curvatures) & tip1_curvatures > 0 & tip1_curvatures < 100;

    if sum(valid_velocity_mask) > 0 || sum(valid_curvature_mask) > 0
        generate_tip_detailed_analysis(valid_times, tip1_velocities, tip1_curvatures, 1);
    else
        fprintf('尖端1数据不足，跳过专项分析\n');
    end
else
    fprintf('尖端1数据不足，跳过专项分析\n');
end

%--- 尖端2专项分析
fprintf('\n进行尖端2专项分析...\n');
if length(valid_times) >= 2 && size(tip_velocity_history_6, 2) >= 2 && size(tip_curvature_history_6, 2) >= 2
    tip2_velocities = tip_velocity_history_6(valid_time_steps, 2);
    tip2_curvatures = tip_curvature_history_6(valid_time_steps, 2);

    % 数据验证
    valid_velocity_mask = abs(tip2_velocities) > 1e-6;
    valid_curvature_mask = isfinite(tip2_curvatures) & tip2_curvatures > 0 & tip2_curvatures < 100;

    if sum(valid_velocity_mask) > 0 || sum(valid_curvature_mask) > 0
        generate_tip_detailed_analysis(valid_times, tip2_velocities, tip2_curvatures, 2);
    else
        fprintf('尖端2数据不足，跳过专项分析\n');
    end
else
    fprintf('尖端2数据不足，跳过专项分析\n');
end

%--- 尖端3专项分析
fprintf('\n进行尖端3专项分析...\n');
if length(valid_times) >= 2 && size(tip_velocity_history_6, 2) >= 3 && size(tip_curvature_history_6, 2) >= 3
    tip3_velocities = tip_velocity_history_6(valid_time_steps, 3);
    tip3_curvatures = tip_curvature_history_6(valid_time_steps, 3);

    % 数据验证
    valid_velocity_mask = abs(tip3_velocities) > 1e-6;
    valid_curvature_mask = isfinite(tip3_curvatures) & tip3_curvatures > 0 & tip3_curvatures < 100;

    if sum(valid_velocity_mask) > 0 || sum(valid_curvature_mask) > 0
        generate_tip_detailed_analysis(valid_times, tip3_velocities, tip3_curvatures, 3);
    else
        fprintf('尖端3数据不足，跳过专项分析\n');
    end
else
    fprintf('尖端3数据不足，跳过专项分析\n');
end

%--- 尖端5专项分析
fprintf('\n进行尖端5专项分析...\n');
if length(valid_times) >= 2 && size(tip_velocity_history_6, 2) >= 5 && size(tip_curvature_history_6, 2) >= 5
    tip5_velocities = tip_velocity_history_6(valid_time_steps, 5);
    tip5_curvatures = tip_curvature_history_6(valid_time_steps, 5);

    % 数据验证
    valid_velocity_mask = abs(tip5_velocities) > 1e-6;
    valid_curvature_mask = isfinite(tip5_curvatures) & tip5_curvatures > 0 & tip5_curvatures < 100;

    if sum(valid_velocity_mask) > 0 || sum(valid_curvature_mask) > 0
        generate_tip_detailed_analysis(valid_times, tip5_velocities, tip5_curvatures, 5);
    else
        fprintf('尖端5数据不足，跳过专项分析\n');
    end
else
    fprintf('尖端5数据不足，跳过专项分析\n');
end

%--- 尖端6专项分析
fprintf('\n进行尖端6专项分析...\n');
if length(valid_times) >= 2 && size(tip_velocity_history_6, 2) >= 6 && size(tip_curvature_history_6, 2) >= 6
    tip6_velocities = tip_velocity_history_6(valid_time_steps, 6);
    tip6_curvatures = tip_curvature_history_6(valid_time_steps, 6);

    % 数据验证
    valid_velocity_mask = abs(tip6_velocities) > 1e-6;
    valid_curvature_mask = isfinite(tip6_curvatures) & tip6_curvatures > 0 & tip6_curvatures < 100;

    if sum(valid_velocity_mask) > 0 || sum(valid_curvature_mask) > 0
        generate_tip_detailed_analysis(valid_times, tip6_velocities, tip6_curvatures, 6);
    else
        fprintf('尖端6数据不足，跳过专项分析\n');
    end
else
    fprintf('尖端6数据不足，跳过专项分析\n');
end

%--- 保存最终数据
save_filename = sprintf('time_controlled_phase_field_%.1fs_%s.mat', total_time, datestr(now, 'yyyymmdd_HHMMSS'));
save(save_filename, ...
     'phi_history', 'U_history', 'time_history', ...
     'interface_history', 'contour_history', ...
     'tip_position_history_6', 'tip_distance_history_6', 'tip_velocity_history_6', 'tip_curvature_history_6', ...
     'six_tip_angles', 'six_tip_directions', 'six_reference_lines', ...
     'tip_history_buffer_6', ...
     'phi', 'U', 'Nx', 'Ny', 'dx', 'dy', 'total_time', 'dtime', 'current_time', 'save_count', ...
     'print_interval', 'max_steps', 'k_partition', 'epsilon_aniso', 'm_aniso', 'omega0_aniso', ...
     'lambda_coup', 'D_coefficient', 'theta_field');

%--- 计算总运行时间
compute_time = etime(clock(), time0);
fprintf('总计算时间: %.2f 秒\n', compute_time);

function [ip, im, jp, jm] = apply_boundary_conditions(i, j, Nx, Ny, bc_type)
    % 应用边界条件：周期性或零梯度
    % 支持两种边界条件类型以提供更大的灵活性

    jp = j + 1;    % 上方邻居点
    jm = j - 1;    % 下方邻居点
    ip = i + 1;    % 右方邻居点
    im = i - 1;    % 左方邻居点

    switch lower(bc_type)
        case 'periodic'
            % 周期性边界条件（原有逻辑）
            % 左右边界周期连接
            if im == 0,        im = Nx;    end    % 左边界映射到右边界
            if ip == (Nx+1),   ip = 1;     end    % 右边界映射到左边界
            % 上下边界周期连接
            if jm == 0,        jm = Ny;    end    % 下边界映射到上边界
            if jp == (Ny+1),   jp = 1;     end    % 上边界映射到下边界

        case 'zero_gradient'
            % 零梯度边界条件（新增逻辑）
            % 边界处使用最近的内点值，实现梯度为零的条件
            % 左右边界零梯度
            if im == 0,        im = 1;     end    % 左边界使用第一个内点
            if ip == (Nx+1),   ip = Nx;    end    % 右边界使用最后一个内点
            % 上下边界零梯度
            if jm == 0,        jm = 1;     end    % 下边界使用第一个内点
            if jp == (Ny+1),   jp = Ny;    end    % 上边界使用最后一个内点

        otherwise
            error('未知的边界条件类型: %s。请使用 ''periodic'' 或 ''zero_gradient''', bc_type);
    end
end

%% ================== 可视化函数 ==================

function initialize_corrected_visualization(Nx, Ny)
    figure('Name', '修正相场-浓度场耦合模拟', 'Position', [100, 100, 1800, 500]);

    % 相场显示
    subplot(1, 3, 1);
    imagesc(zeros(Ny, Nx));
    colormap(jet);
    colorbar;
    title('相场分布 ϕ');
    xlabel('x'); ylabel('y');
    axis equal tight;
    caxis([-1, 1]);  % 固定相场显示范围

    % 浓度场显示
    subplot(1, 3, 2);
    imagesc(zeros(Ny, Nx));
    colormap(jet);
    colorbar;
    title('浓度场分布 U');
    xlabel('x'); ylabel('y');
    axis equal tight;
    caxis([0, 2.5]);  % 固定浓度场显示范围

    % C/C∞云图显示
    subplot(1, 3, 3);
    imagesc(zeros(Ny, Nx));
    colormap(jet);
    colorbar;
    title('浓度比分布 C/C∞');
    xlabel('x'); ylabel('y');
    axis equal tight;
    caxis([0, 2.5]);  % 固定C/C∞显示范围

    drawnow;
end

function update_corrected_visualization(phi, U, istep, probe_data, Nx, Ny, y_ratio, k_partition)
    % 更新可视化，集成探测线分析功能
    %
    % 输入:
    %   phi, U - 相场和浓度场
    %   istep - 当前步数
    %   probe_data - 探测线数据(可选)
    %   Nx, Ny - 网格尺寸
    %   y_ratio - 探测线位置比例
    %   k_partition - 溶质分配系数

    figure(1);

    % 更新相场显示
    subplot(1, 3, 1);
    imagesc(phi');
    title(sprintf('相场分布 ϕ - 步数 %d', istep));
    caxis([-1, 1]);  % 保持固定范围以捕捉细微变化
    xlabel('x'); ylabel('y');
    axis equal tight;
    colorbar;  % 确保颜色条图例存在

    % 更新浓度场显示
    subplot(1, 3, 2);
    imagesc(U');
    title(sprintf('浓度场分布 U - 步数 %d', istep));
    caxis([0, 2.5]);  % 保持固定范围以捕捉细微变化
    xlabel('x'); ylabel('y');
    axis equal tight;
    colorbar;  % 确保颜色条图例存在

    % 计算并更新C/C∞云图显示
    subplot(1, 3, 3);
    CCinf_field = calculate_CCinf_from_fields(phi, U, k_partition);
    imagesc(CCinf_field');
    title(sprintf('浓度比分布 C/C∞ - 步数 %d', istep));
    caxis([0, 2.5]);  % 保持固定范围以捕捉细微变化
    xlabel('x'); ylabel('y');
    axis equal tight;
    colorbar;  % 确保颜色条图例存在

    drawnow;
end

%% ================== 数据探测线分析函数 ==================

function CCinf_ratio = calculate_CCinf_from_fields(phi, U, k_partition)
    % 计算C/C∞比值
    % 基于公式: U = ((2C/C∞)/(1+k-(1-k)φ) - 1)/(1-k)
    % 反解得到: C/C∞ = (U*(1-k) + 1) * (1+k-(1-k)φ) / 2
    %
    % 输入:
    %   phi - 相场值
    %   U - 浓度场值
    %   k_partition - 溶质分配系数
    %
    % 输出:
    %   CCinf_ratio - C/C∞比值

    CCinf_ratio = (U .* (1 - k_partition) + 1) .* ...
                  (1 + k_partition - (1 - k_partition) .* phi) / 2;

    % 添加数值稳定性处理
    CCinf_ratio(isnan(CCinf_ratio)) = 0;
    CCinf_ratio(isinf(CCinf_ratio)) = 0;
end

function [y_indices, CCinf_data, phi_line, U_line] = extract_probe_data(phi, U, k_partition, Nx, Ny, x_ratio)
    % 在指定x位置提取垂直探测线数据
    %
    % 输入:
    %   phi, U - 相场和浓度场
    %   k_partition - 溶质分配系数
    %   Nx, Ny - 网格尺寸
    %   x_ratio - x位置比例(0到1之间)
    %
    % 输出:
    %   y_indices - y网格序号
    %   CCinf_data - C/C∞比值数据
    %   phi_line - 探测线上的相场值
    %   U_line - 探测线上的浓度场值

    % 计算探测线的x坐标
    x_probe = round(Nx * x_ratio);
    if x_probe < 1, x_probe = 1; end
    if x_probe > Nx, x_probe = Nx; end

    % 提取垂直探测线上的数据（固定x，遍历y）
    phi_line = phi(:, x_probe);
    U_line = U(:, x_probe);

    % 计算C/C∞比值
    CCinf_data = calculate_CCinf_from_fields(phi_line, U_line, k_partition);

    % y坐标为网格序号
    y_indices = 1:Ny;
end

function plot_probe_analysis_overlay(phi, U, CCinf_data, istep, Nx, Ny, y_ratio, k_partition)
    % 在现有可视化基础上叠加探测线分析
    %
    % 输入:
    %   phi, U - 相场和浓度场
    %   CCinf_data - C/C∞比值数据
    %   istep - 当前步数
    %   Nx, Ny - 网格尺寸
    %   y_ratio - 探测线y位置比例
    %   k_partition - 溶质分配系数

    % 获取当前图形窗口
    figure(1);

    % 计算探测线位置
    y_probe = round(Ny * y_ratio);
    x_indices = 1:Nx;

    % 在相场图上添加探测线
    subplot(1, 3, 1);
    hold on;
    plot(x_indices, y_probe * ones(size(x_indices)), 'r-', 'LineWidth', 2);
    title(sprintf('相场分布 ϕ - 步数 %d (探测线y=%d)', istep, y_probe));
    hold off;

    % 在浓度场图上添加探测线
    subplot(1, 3, 2);
    hold on;
    plot(x_indices, y_probe * ones(size(x_indices)), 'r-', 'LineWidth', 2);
    title(sprintf('浓度场分布 U - 步数 %d (探测线y=%d)', istep, y_probe));
    hold off;

    % 在C/C∞云图上添加探测线
    subplot(1, 3, 3);
    hold on;
    plot(x_indices, y_probe * ones(size(x_indices)), 'r-', 'LineWidth', 2);
    title(sprintf('浓度比分布 C/C∞ - 步数 %d (探测线y=%d)', istep, y_probe));
    hold off;

    % 在图形窗口上方添加C/C∞曲线
    % 创建新的图形窗口避免覆盖现有显示
    figure('Name', '探测线分析', 'Position', [300, 100, 800, 400]);

    plot(x_indices, CCinf_data, 'b-', 'LineWidth', 2);
    grid on;
    xlabel('网格序号 (无单位)');
    ylabel('C/C∞ (无单位)');
    title(sprintf('浓度分布探测线分析 - 步数 %d (y=%d)', istep, y_probe));

    % 添加参考线
    yline(1.0, 'r--', 'C∞', 'LineWidth', 1);
    yline(k_partition, 'g--', 'k·C∞', 'LineWidth', 1);

    % 设置固定的y轴范围为0到2.5
    ylim([0, 2.5]);

    drawnow;
end

function save_probe_analysis_data(y_indices, CCinf_data, phi_line, U_line, istep, y_probe, W_over_d0)
    % 保存探测线分析数据
    %
    % 输入:
    %   y_indices - y网格序号
    %   CCinf_data - C/C∞比值数据
    %   phi_line - 探测线相场值
    %   U_line - 探测线浓度场值
    %   istep - 当前步数
    %   y_probe - 探测线y位置

    % 保存为MAT文件
    mat_filename = sprintf('probe_analysis_step_%d_y%d_W_over_d0=%.4f.mat', istep, y_probe, W_over_d0);
    save(mat_filename, 'y_indices', 'CCinf_data', 'phi_line', 'U_line', 'y_probe');

    % 保存为CSV文件
    csv_filename = sprintf('probe_analysis_step_%d_y%d_W_over_d0=%.4f.csv', istep, y_probe, W_over_d0);

    fid = fopen(csv_filename, 'w');
    fprintf(fid, 'x_index,CCinf_ratio,phi,U,y_position\n');
    for i = 1:length(y_indices)
        fprintf(fid, '%d,%.6e,%.6f,%.6f,%d\n', y_indices(i), CCinf_data(i), phi_line(i), U_line(i), y_probe);
    end
    fclose(fid);

    fprintf('探测线数据已保存: %s 和 %s\n', mat_filename, csv_filename);
end

%% ================== 界面追踪可视化函数 ==================

function initialize_interface_visualization(Nx, Ny, dx, dy)
    % Initialize solid-liquid interface tracking visualization window (figure(2))
    % Create 2x2 subplot layout for interface tracking analysis
    %
    % Input:
    %   Nx, Ny - Grid dimensions
    %   dx, dy - Spatial step sizes

    figure(2);
    set(gcf, 'Name', 'Solid-Liquid Interface Tracking Analysis', 'Position', [200, 200, 1200, 900]);

    % Subplot 1: Interface contour real-time display (top-left)
    subplot(2, 2, 1);
    imagesc(zeros(Ny, Nx));
    colormap('gray');
    title('Interface Contour Real-time Display', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('y (grid units)');  % Corrected: because of transpose, x-axis shows y
    ylabel('x (grid units)');  % Corrected: because of transpose, y-axis shows x
    axis equal tight;
    hold on;

    % Initialize interface contour line object
    interface_line = plot([], 'r-', 'LineWidth', 2, 'DisplayName', 'Current Interface');
    hold off;

    % Subplot 2: Interface length evolution curve (top-right)
    subplot(2, 2, 2);
    interface_length_plot = plot([], 'b-', 'LineWidth', 2, 'Marker', 'o');
    title('Interface Length Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time Step');
    ylabel('Interface Length (grid units)');
    grid on;
    xlim([0, 100]);  % Initial range, will be dynamically adjusted
    ylim([0, 1000]); % Initial range, will be dynamically adjusted

    % Subplot 3: Interface geometric feature analysis (bottom-left)
    subplot(2, 2, 3);
    scatter_plot = scatter([], [], 20, 'filled');
    title('Interface Coordinate Distribution', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('x Coordinate (grid units)');
    ylabel('y Coordinate (grid units)');
    grid on;
    axis equal;
    xlim([0, Nx*dx]);
    ylim([0, Ny*dy]);

    % Subplot 4: Statistical information panel (bottom-right)
    subplot(2, 2, 4);
    axis off;
    info_text = text(0.1, 0.9, '', 'FontSize', 10, 'VerticalAlignment', 'top');
    title('Statistical Information', 'FontSize', 12, 'FontWeight', 'bold');

    % Add main title
    sgtitle('Solid-Liquid Interface Tracking Dynamic Analysis System', 'FontSize', 14, 'FontWeight', 'bold');

    drawnow;
end

function update_interface_visualization(phi, contour_data, interface_length, istep, time_history, interface_history, Nx, Ny, dx, dy, tip_pos, tip_distance, tip_velocity, tip_curvature_radius)
    % Update solid-liquid interface tracking visualization window
    %
    % Input parameters:
    %   phi - Current phase field distribution
    %   contour_data - Current interface coordinate data [x_coords, y_coords]
    %   interface_length - Current interface length
    %   istep - Current time step
    %   time_history - Time history array
    %   interface_history - Interface length history array
    %   Nx, Ny - Grid dimensions
    %   dx, dy - Spatial step sizes
    %   tip_pos - Current tip position [x, y]
    %   tip_distance - Distance from tip to center
    %   tip_velocity - Tip growth velocity
    %   tip_curvature_radius - Tip curvature radius

    figure(2);

    % === Subplot 1: Interface contour real-time display ===
    subplot(2, 2, 1);
    imagesc(phi');
    colormap('gray');
    title(sprintf('Interface Contour Real-time Display - Step %d', istep), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('y (grid units)');  % Corrected: because of transpose, x-axis shows y coordinate
    ylabel('x (grid units)');  % Corrected: because of transpose, y-axis shows x coordinate
    axis equal tight;
    hold on;

    % Plot current interface contour - Correct coordinate correspondence
    if ~isempty(contour_data)
        % Because imagesc(phi') is transposed, need to swap x,y coordinates of contour_data
        plot(contour_data(:, 2)/dx, contour_data(:, 1)/dy, 'r-', 'LineWidth', 2, 'DisplayName', 'Current Interface (φ=0)');
    end
    hold off;

    % === Subplot 2: Interface length evolution curve ===
    subplot(2, 2, 2);
    if length(interface_history) > 1
        valid_indices = interface_history > 0;
        plot((0:length(interface_history)-1)*(50), interface_history, 'b-', 'LineWidth', 2, 'Marker', 'o');
    else
        plot([0], [interface_length], 'b-', 'LineWidth', 2, 'Marker', 'o');
    end

    title('Interface Length Evolution', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time Step');
    ylabel('Interface Length (grid units)');
    grid on;

    % Add current value annotation
    text(0.02, 0.98, sprintf('Current: %.1f', interface_length), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'FontSize', 10);

    % === Subplot 3: Interface contour scatter plot ===
    subplot(2, 2, 3);
    if ~isempty(contour_data)
        scatter(contour_data(:, 2), contour_data(:, 1), 15, contour_data(:, 1), 'filled');
        colormap('jet');
        colorbar;
        ylabel(colorbar, 'x coordinate (grid units)');  % Add colorbar label to explain color meaning
    end
    title('Interface Coordinate Point Distribution', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('y coordinate (grid units)');  % Corrected: consistent with subplot 1
    ylabel('x coordinate (grid units)');  % Corrected: consistent with subplot 1
    grid on;
    axis equal;
    xlim([0, Ny*dy]);  % Corrected: coordinate range after transpose
    ylim([0, Nx*dx]);  % Corrected: coordinate range after transpose

    % === Subplot 4: Statistical information panel ===
    subplot(2, 2, 4);
    cla;         % Clear subplot content but keep settings
    hold off;    % Ensure no hold state
    axis off;    % Reset axis off for clean display environment

    % Calculate statistical information
    interface_points = size(contour_data, 1);
    solid_fraction = sum(phi(:) > 0) / (Nx * Ny) * 100;
    liquid_fraction = 100 - solid_fraction;

    % Estimate growth rate (if historical data available)
    if length(interface_history) > 1
        recent_growth = interface_history(end) - interface_history(max(1, end-5));
        growth_rate = recent_growth / min(5, length(interface_history)-1) * 50; % Growth rate per 50 steps
    else
        growth_rate = 0;
    end

    % Display statistical information - Check first then build to avoid array index errors
    if isempty(tip_pos) || length(tip_pos) < 2
        % If no tip detected or data incomplete, display N/A information
        info_str = {
            sprintf('========== Interface Tracking Statistics ==========');
            sprintf('Current Time Step: %d', istep);
            sprintf('Total Interface Length: %.2f grid units', interface_length);
            sprintf('Interface Points: %d', interface_points);
            sprintf('');
            sprintf('========== Phase Distribution Statistics ==========');
            sprintf('Solid Fraction: %.2f%%', solid_fraction);
            sprintf('Liquid Fraction: %.2f%%', liquid_fraction);
            sprintf('');
            sprintf('========== Tip Dynamics Information ==========');
            sprintf('Tip Position: N/A');
            sprintf('Tip Distance from Center: N/A');
            sprintf('Tip Velocity: N/A');
            sprintf('Tip Curvature Radius: N/A');
            sprintf('');
            sprintf('========== Kinetic Information ==========');
            sprintf('Growth Rate: %.3f step⁻¹', growth_rate);
            sprintf('Interface Threshold: φ = 0 (adapted for [-1,1] range)');
            sprintf('');
            sprintf('========== Grid Information ==========');
            sprintf('Grid Size: %d × %d', Nx, Ny);
            sprintf('Spatial Step: %.2f × %.2f', dx, dy);
            sprintf('Computational Domain: %.1f × %.1f', Nx*dx, Ny*dy);
        };
    else
        % Display when valid tip data available
        curvature_str = 'N/A';
        if ~isempty(tip_curvature_radius) && tip_curvature_radius < 1000
            curvature_str = sprintf('%.2f grid units', tip_curvature_radius);
        end

        info_str = {
            sprintf('========== Interface Tracking Statistics ==========');
            sprintf('Current Time Step: %d', istep);
            sprintf('Total Interface Length: %.2f grid units', interface_length);
            sprintf('Interface Points: %d', interface_points);
            sprintf('');
            sprintf('========== Phase Distribution Statistics ==========');
            sprintf('Solid Fraction: %.2f%%', solid_fraction);
            sprintf('Liquid Fraction: %.2f%%', liquid_fraction);
            sprintf('');
            sprintf('========== Tip Dynamics Information ==========');
            sprintf('Tip Position: (%.1f, %.1f)', tip_pos(1), tip_pos(2));
            sprintf('Tip Distance from Center: %.2f grid units', tip_distance);
            sprintf('Tip Velocity: %.4f grid units/s', tip_velocity);
            sprintf('Tip Curvature Radius: %s', curvature_str);
            sprintf('');
            sprintf('========== Kinetic Information ==========');
            sprintf('Growth Rate: %.3f step⁻¹', growth_rate);
            sprintf('Interface Threshold: φ = 0 (adapted for [-1,1] range)');
            sprintf('');
            sprintf('========== Grid Information ==========');
            sprintf('Grid Size: %d × %d', Nx, Ny);
            sprintf('Spatial Step: %.2f × %.2f', dx, dy);
            sprintf('Computational Domain: %.1f × %.1f', Nx*dx, Ny*dy);
        };
    end

    % Display statistical text with proper centering
    text(0.5, 0.5, info_str, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8);

    % Display title
    title('Statistical Information', 'FontSize', 12, 'FontWeight', 'bold');

    drawnow;
end

% ================== 视频管理函数 ==================

function video_handle = video_manager(action, filename, fig_handle)
    % Video management function for recording simulation animations
    %
    % Inputs:
    %   action    - 'init', 'add', or 'close'
    %   filename  - video filename (for 'init' action)
    %   fig_handle - figure handle to record (for 'add' action)
    %
    % Output:
    %   video_handle - handle to the video writer object

    persistent video_writer;

    switch action
        case 'init'
            % Initialize video writer
            video_writer = VideoWriter(filename, 'MPEG-4');
            video_writer.FrameRate = 10;
            video_writer.Quality = 90;
            open(video_writer);
            video_handle = video_writer;

        case 'add'
            % Add current frame to video
            if isempty(video_writer)
                error('Video writer not initialized. Call video_manager(''init'', ...) first.');
            end
            frame = getframe(fig_handle);
            writeVideo(video_writer, frame);
            video_handle = video_writer;

        case 'close'
            % Close video writer
            if ~isempty(video_writer)
                close(video_writer);
                video_writer = [];
            end
            video_handle = [];

        otherwise
            error('Unknown action: %s. Valid actions are ''init'', ''add'', ''close''.', action);
    end
end
