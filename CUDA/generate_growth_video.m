%% 生成动态生长视频（φ 和 C/C∞）
%
% 用法：
%   generate_growth_video('frames/', 1000, 1000);
%
% 功能：
%   读取所有帧文件，生成动态生长视频
%   显示相场 φ 和归一化浓度 C/C∞ 的演化

function generate_growth_video(frames_dir, Nx, Ny)
    fprintf('========================================\n');
    fprintf('生成动态生长视频\n');
    fprintf('========================================\n\n');

    %% 1. 查找所有帧文件
    fprintf('[步骤 1/5] 扫描帧文件...\n');

    phi_files = dir(fullfile(frames_dir, 'phi_frame_*.bin'));
    U_files = dir(fullfile(frames_dir, 'U_frame_*.bin'));

    if isempty(phi_files) || isempty(U_files)
        error('未找到帧文件！请先运行CUDA模拟生成中间结果。');
    end

    num_frames = length(phi_files);
    fprintf('找到 %d 帧数据\n', num_frames);

    % 按文件名排序（确保顺序正确）
    [~, idx] = sort({phi_files.name});
    phi_files = phi_files(idx);
    [~, idx] = sort({U_files.name});
    U_files = U_files(idx);

    %% 2. 读取所有帧数据到内存
    fprintf('\n[步骤 2/5] 读取帧数据...\n');

    % 预分配数组
    phi_history = zeros(Ny, Nx, num_frames, 'single');
    U_history = zeros(Ny, Nx, num_frames, 'single');

    % 读取每一帧
    for i = 1:num_frames
        % 读取相场
        fid = fopen(fullfile(frames_dir, phi_files(i).name), 'r');
        temp = fread(fid, [Nx, Ny], 'float32');
        phi_history(:, :, i) = temp';
        fclose(fid);

        % 读取浓度场
        fid = fopen(fullfile(frames_dir, U_files(i).name), 'r');
        temp = fread(fid, [Nx, Ny], 'float32');
        U_history(:, :, i) = temp';
        fclose(fid);

        if mod(i, 10) == 0 || i == num_frames
            fprintf('  已读取 %d/%d 帧...\n', i, num_frames);
        end
    end

    fprintf('所有帧数据读取完成！\n');

    %% 3. 计算 C/C∞
    fprintf('\n[步骤 3/5] 计算 C/C∞...\n');

    k_partition = 0.15;
    CCinf_history = ((1 + k_partition - (1 - k_partition) .* phi_history) / 2) .* (1 + U_history);
    CCinf_history(isnan(CCinf_history)) = 0;
    CCinf_history(isinf(CCinf_history)) = 0;

    fprintf('C/C∞ 计算完成！\n');

    %% 4. 创建视频
    fprintf('\n[步骤 4/5] 生成视频...\n');

    % 创建视频写入器
    video_phi = VideoWriter('dendritic_growth_phi.mp4', 'MPEG-4');
    video_CCinf = VideoWriter('dendritic_growth_CCinf.mp4', 'MPEG-4');

    video_phi.FrameRate = 10;  % 每秒10帧
    video_CCinf.FrameRate = 10;

    open(video_phi);
    open(video_CCinf);

    % 生成每一帧
    figure('Position', [50, 50, 1800, 800], 'Visible', 'off');

    for i = 1:num_frames
        % 提取当前帧
        phi_frame = phi_history(:, :, i);
        U_frame = U_history(:, :, i);
        CCinf_frame = CCinf_history(:, :, i);

        % ========== 创建复合图像 ==========
        clf;

        % 子图1: 相场 φ
        subplot(1, 2, 1);
        imagesc(phi_frame);
        colormap(gca, 'jet');
        caxis([-1, 1]);
        colorbar;
        title(sprintf('相场分布 φ (帧 %d/%d)', i, num_frames), 'Interpreter', 'none');
        xlabel('x'); ylabel('y');
        axis equal tight;

        % 子图2: C/C∞
        subplot(1, 2, 2);
        imagesc(CCinf_frame);
        colormap(gca, 'jet');
        caxis([0, 2.5]);  % 固定颜色范围
        colorbar;
        title(sprintf('归一化浓度 C/C∞ (帧 %d/%d)', i, num_frames), 'Interpreter', 'none');
        xlabel('x'); ylabel('y');
        axis equal tight;

        % 添加时间信息
        annotation('textbox', [0.02, 0.92, 0.15, 0.05], 'String', ...
                  sprintf('时间: %.2fτ₀', (i-1) * 0.5), ...
                  'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', ...
                  'BackgroundColor', 'white', 'Units', 'normalized');

        % 保存当前帧到视频
        frame = getframe(gcf);
        writeVideo(video_phi, frame);
        writeVideo(video_CCinf, frame);

        if mod(i, 20) == 0 || i == num_frames
            fprintf('  已处理 %d/%d 帧...\n', i, num_frames);
        end
    end

    close(video_phi);
    close(video_CCinf);

    fprintf('视频生成完成！\n');
    fprintf('  相场视频: dendritic_growth_phi.mp4\n');
    fprintf('  C/C∞ 视频: dendritic_growth_CCinf.mp4\n');

    %% 5. 生成综合对比图（所有帧的缩略图）
    fprintf('\n[步骤 5/5] 生成时间演化缩略图...\n');

    create_timeline_summary(phi_history, CCinf_history, num_frames, Nx, Ny);

    fprintf('\n========================================\n');
    fprintf('视频生成完成！\n');
    fprintf('========================================\n');
end

%% 创建时间演化摘要图
function create_timeline_summary(phi_history, CCinf_history, num_frames, Nx, Ny)
    fprintf('生成时间演化摘要图...\n');

    figure('Position', [100, 100, 2000, 600], 'Name', 'Timeline Summary');

    % 选择关键帧显示
    num_key_frames = min(10, num_frames);
    key_indices = round(linspace(1, num_frames, num_key_frames));

    % 相场时间线
    for i = 1:num_key_frames
        idx = key_indices(i);
        subplot(2, num_key_frames, i);
        imagesc(phi_history(:, :, idx));
        colormap(gca, 'jet');
        caxis([-1, 1]);
        title(sprintf('%.1fτ₀', (idx-1) * 0.5), 'Interpreter', 'none');
        axis equal tight;
        axis off;
    end

    % C/C∞ 时间线
    for i = 1:num_key_frames
        idx = key_indices(i);
        subplot(2, num_key_frames, num_key_frames + i);
        imagesc(CCinf_history(:, :, idx));
        colormap(gca, 'turbo');
        caxis([0, 1.5]);
        title(sprintf('%.1fτ₀', (idx-1) * 0.5), 'Interpreter', 'none');
        axis equal tight;
        axis off;
    end

    sgtitle('相场与 C/C∞ 时间演化', 'FontSize', 14, 'FontWeight', 'bold');

    % 保存
    saveas(gcf, 'timeline_summary.png');
    fprintf('  时间演化摘要图已保存: timeline_summary.png\n');
end
