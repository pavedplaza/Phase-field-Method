%% 可视化模拟结果（三场：φ, U, C/C∞）
%
% 用法：view_results
%
% 输出：三张2D分布图（相场、浓度场、归一化浓度）

function view_results()
    Nx = 1000;
    Ny = 1000;
    build_dir = 'build';

    % 读取数据
    fprintf('读取数据...\n');
    fid = fopen(fullfile(build_dir, 'phi_final.bin'), 'r');
    phi = fread(fid, [Nx, Ny], 'float32')';
    fclose(fid);

    fid = fopen(fullfile(build_dir, 'U_final.bin'), 'r');
    U = fread(fid, [Nx, Ny], 'float32')';
    fclose(fid);

    % 计算 C/C∞
    k_partition = 0.15;
    CCinf = ((1 + k_partition - (1 - k_partition) .* phi) / 2) .* (1 + U);

    % 保存到工作区
    assignin('base', 'phi', phi);
    assignin('base', 'U', U);
    assignin('base', 'CCinf', CCinf);
    fprintf('数据已保存到工作区\n');

    % 绘图（3个子图）
    figure('Position', [50, 50, 1800, 500], 'Name', 'Phase Field Results');

    % 子图1: 相场 φ
    subplot(1, 3, 1);
    imagesc(phi);
    colormap(gca, 'jet');
    caxis([-1, 1]);
    colorbar;
    title('相场 φ', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('x'); ylabel('y');
    axis equal tight;

    % 子图2: 浓度场 U
    subplot(1, 3, 2);
    imagesc(U);
    colormap(gca, 'jet');
    colorbar;
    title('浓度场 U', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('x'); ylabel('y');
    axis equal tight;

    % 子图3: 归一化浓度 C/C∞
    subplot(1, 3, 3);
    imagesc(CCinf);
    colormap(gca, 'jet');
    caxis([0, 2.5]);
    colorbar;
    title('归一化浓度 C/C∞', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('x'); ylabel('y');
    axis equal tight;

    % 保存
    saveas(gcf, 'results.png');
    fprintf('\n图像已保存: results.png\n');
end
