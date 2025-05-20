clc;
clear;

% 自定义 x 轴标签
% 定义要标注的位置和对应的标签
positions = [1, 5, 9, 12, 14, 17, 21, 25];
labels = {'-3M', '-2M', '-1M', '-8D', '8D', '1M', '2M', '3M'};

% 设置默认字体为 Times New Roman
set(groot, 'defaultAxesFontName', 'Times New Roman');

% 指定路径
path = 'E:\Data\站点数据_1217\绘图数据';

% 获取路径下的所有文件
files = dir(fullfile(path, '*.*')); % *.* 表示匹配所有文件
files = files(~[files.isdir]); % 排除子文件夹

% 遍历每个文件
for f = 2:length(files)
    % 获取完整文件路径
    file_name_first = fullfile(files(f).folder, files(f).name);
    load(file_name_first);
    table_first = combined_result_table;
    
    % 第二个文件路径和加载
    file_name_other = fullfile(files(4).folder, files(4).name);
    load(file_name_other);
    table_other = combined_result_table;
    
    % 拼接两个文件的表数据
    combined_table = [table_first; table_other]; % 垂直拼接

%     GPP_f_data = combined_table.('GPP_DT''uWUE_T');
%     GPP_f_data = GPP_f_data(:)'; 

    column_name = 'GPP_DT';
    table_data = combined_table.(column_name);
    table_data = table_data(:)';
    
    % 指定中心点和范围
    center_idx = 40;
    range = 12;

    % 计算索引范围
    start_idx = center_idx - range;
    end_idx = center_idx + range;

    % 提取数据范围
    selected_data = table_data(start_idx:end_idx);

    % 提取文件名部分（不包括路径和扩展名）
    [~, base_name, ~] = fileparts(file_name_first);
    words = split(base_name, '_');
    % 获取第三个和第四个单词
    year_word = words{4};
    site_word = words{5};

    % 构建新的输出文件名
    output_file = sprintf('%s_%s_%s_new.svg', year_word, site_word, column_name);

    % 根据 fourth_word 是否包含特定关键字选择颜色
    if contains(column_name, 'ESI', 'IgnoreCase', true)
        color_fast = '#F12C32';
        color_mean = '#FB6C3F';
        color_slow = '#FBA250'; 
    elseif contains(column_name, 'LAI', 'IgnoreCase', true)
        color_fast = '#27A09E';
        color_mean = '#30CE88';
        color_slow = '#7DE393';  
    elseif contains(column_name, 'SIF', 'IgnoreCase', true)
        color_fast = '#4FFD95';
        color_mean = '#97FD52';
        color_slow = '#C5FB2E';
    elseif contains(column_name, '_T', 'IgnoreCase', true)
        color_mean = [12/255, 114/255, 31/255];
    elseif contains(column_name, 'NEP', 'IgnoreCase', true)
        color_fast = '#11AB0B';
        color_mean = '#17D80F';
        color_slow = '#72F46D';
    elseif contains(column_name, 'GPP', 'IgnoreCase', true)
        color_mean = [9/255, 91/255, 208/255];
    elseif contains(column_name, 'NEE', 'IgnoreCase', true)
        color_mean = [23/255, 216/255, 15/255];
    end

    period = 1:25;
    serious_period = 13;

    % 创建新的图窗并设置其位置和大小
    figure('Position', [800, 100, 800, 580]);
    
    % 绘制折线图
%     plot(period, selected_GPP_f_data, '-','Color', color_fast, 'LineWidth', 4);
    plot(period, selected_data, '-','Color', color_mean, 'LineWidth', 4);
    hold on;
    
    period_points = readmatrix("E:\Data\站点数据_1217\千烟洲\FD_178\FD_178_periods.csv");

    % 根据 fourth_word 包含的关键字选择设置
    if contains(column_name, 'ESI')
        xlim_custom = [0 26];
        xticks_custom = 0:1:26;
        ylim_custom = [-1.5 1.5];
        yticks_custom = -1.5:0.5:1.5;
    elseif contains(column_name, 'LAI')
        xlim_custom = [0 26];
        xticks_custom = 0:1:26;
        ylim_custom = [0 6];
        yticks_custom = 0:1:6;
    elseif contains(column_name, 'SIF')
        xlim_custom = [0 26];
        xticks_custom = 0:1:26;
        ylim_custom = [-0.2 1];
        yticks_custom = -0.2:0.2:1;
    elseif contains(column_name, '_T')
        xlim_custom = [0 26];
        xticks_custom = 0:1:26;
        ylim_custom = [-0.5 5.5];
        yticks_custom = -0.5:1:5.5;
    elseif contains(column_name, 'NEP')
        xlim_custom = [0 26];
        xticks_custom = 0:1:26;
        ylim_custom = [-1.5 2];
        yticks_custom = -1.5:0.5:2;
    elseif contains(column_name, 'GPP')
        xlim_custom = [0 26];        
        xticks_custom = 0:1:26;      
        ylim_custom = [-2 10];    
        yticks_custom = -2:2:10;
    elseif contains(column_name, 'NEE')
        xlim_custom = [0 26];
        xticks_custom = 0:1:26;
        ylim_custom = [-6 0];
        yticks_custom = -6:1:0;
    end

    % 设置 x 和 y 轴范围
    xlim(xlim_custom);
    ylim(ylim_custom);
    
    % 设置 x 和 y 轴刻度间隔
    xticks(xticks_custom);
    yticks(yticks_custom);

    if contains(column_name, 'ESI')
        % 设置 y 轴刻度，只显示部分刻度
        set(gca, 'YTick', [-1.0, -0.5, 0, 0.5, 1]);
    elseif contains(column_name, 'LAI')
        set(gca, 'YTick', [1, 2, 3, 4, 5]);
    elseif contains(column_name, 'SIF')
        set(gca, 'YTick', [0, 0.2, 0.4, 0.6, 0.8]);
    elseif contains(column_name, '_T')
        set(gca, 'YTick', [0, 1.0, 2.0, 3.0, 4.0, 5.0]);
    elseif contains(column_name, 'NEP')
        set(gca, 'YTick', [-1.0, -0.5, 0, 0.5, 1.0, 1.5]);
    elseif contains(column_name, 'GPP')
        set(gca, 'YTick', [0, 2, 4, 6, 8]);
    elseif contains(column_name, 'NEE')
        set(gca, 'YTick', [-5, -4, -3, -2, -1]);
    end
    
    % 获取当前的 y 轴刻度
    yTicks = get(gca, 'YTick');
    
    % 将刻度四舍五入到一位小数，并转换成字符串格式
    formattedLabels = sprintf('%.1f\n', yTicks);
    
    % 设置 y 轴刻度的文本格式
    set(gca, 'YTickLabel', formattedLabels);
    
    % 创建阴影区域
    x_fill_shadow_1 = [period_points(1), serious_period, serious_period, period_points(1)];
    y_fill_shadow_1 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
    
    % 添加阴影
    fill(x_fill_shadow_1, y_fill_shadow_1, [128/255, 128/255, 128/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % 创建阴影区域
    x_fill_shadow_2 = [serious_period, period_points(2), period_points(2), serious_period];
    y_fill_shadow_2 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
    
    % 添加阴影
    fill(x_fill_shadow_2, y_fill_shadow_2, [128/255, 128/255, 128/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % 创建阴影区域
    x_fill_shadow_3 = [period_points(2), 26, 26, period_points(2)];
    y_fill_shadow_3 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
    
    % 添加阴影
    fill(x_fill_shadow_3, y_fill_shadow_3, [128/255, 128/255, 128/255], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % 再次绘制
%     plot(period, selected_GPP_f_data, '-','Color', color_fast, 'LineWidth', 4);
    plot(period, selected_data, '-','Color', color_mean, 'LineWidth', 4);

    if contains(column_name, 'ESI')
        title_str = sprintf('%s Response to Flash Drought', 'ESI');
        ylabel_str = '';
    elseif contains(column_name, 'LAI')
        title_str = sprintf('%s Response to Flash Drought', 'LAI');
        ylabel_str = '';
    elseif contains(column_name, 'SIF')
        title_str = sprintf('%s Response to Flash Drought', 'SIF');
        ylabel_str = 'mWm⁻² nm⁻¹ sr⁻¹';
    elseif contains(column_name, '_T')
        title_str = sprintf('%s Response to Flash Drought Based on Site', 'T');
        ylabel_str = 'T (mm day⁻¹)';
    elseif contains(column_name, 'NEP')
        title_str = sprintf('%s Response to Flash Drought', 'NEP');
        ylabel_str = 'g C m⁻² d⁻¹';
    elseif contains(column_name, 'GPP')
        title_str = sprintf('%s Response to Flash Drought Based on Site', 'GPP');
        ylabel_str = 'GPP (g C m⁻² d⁻¹)';
    elseif contains(column_name, 'NEE')
        title_str = sprintf('%s Response to Flash Drought Based on Site', 'NEE');
        ylabel_str = 'g C m⁻² d⁻¹';
    end
 
    % 获取当前坐标轴对象
    ax = gca;
    
    % 设置 x 轴和 y 轴的线条粗细
    ax.XAxis.LineWidth = 2;
    ax.YAxis.LineWidth = 2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    ax.YAxis.FontWeight = 'bold';
    
    
    set(gca,  'XTickLabel', []);

    % 获取坐标轴位置信息
    ax_position = ax.Position;

    % 设置标题
    title_handle_1 = title(title_str, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');

    % 设置标签
    ylabel(ylabel_str, 'FontName', 'Times New Roman', 'FontSize', 22,'FontWeight', 'bold');
    xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');

    if contains(column_name, 'GPP')
        legend_handle = legend('CN-Qia Site', 'Location', 'northwest','FontSize', 22, 'FontWeight', 'bold');
    elseif contains(column_name, 'NEE')
        legend_handle = legend('CN-Qia Site', 'Location', 'northwest','FontSize', 22, 'FontWeight', 'bold');
    elseif contains(column_name, '_T')
        legend_handle = legend('CN-Qia Site', 'Location', 'northwest','FontSize', 22, 'FontWeight', 'bold');
    end

    legend_handle.Box = 'off';

    ax.YAxis.Label.Position(1) = -2;

    % 根据不同的 fourth_word 来调整添加额外的垂直线条和调整 x 轴标签位置
    switch column_name
        case 'ESI'

            y_norm2 = 0.05 / (ylim_custom(2) - ylim_custom(1));
            y_relative = y_norm2 * ax_position(4);
            abs_length = y_relative*((ylim_custom(2)-ylim_custom(1))/ax_position(4));

            % 在原有刻度上添加额外的标记
            for i = 1:length(positions)
                x_position = positions(i);
                x_norm2 = (x_position - xlim_custom(1)) / (xlim_custom(2) - xlim_custom(1));

                x_relative = x_norm2 * ax_position(3);

                % 绘制垂直于 x 轴的线条
                line([x_relative*(26/ax_position(3)), x_relative*(26/ax_position(3))], [abs_length + ylim_custom(1),ylim_custom(1)], ...
                    'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            
            % 调整 x 轴标签位置
            y_norm3 = 0.15 / (ylim_custom(2) - ylim_custom(1));
            y_relative3 = y_norm3 * ax_position(4);
            abs_length3 = y_relative3*((ylim_custom(2)-ylim_custom(1))/ax_position(4));

            ax.XAxis.Label.Position(2) = ax.XAxis.Label.Position(2) - abs_length3;
        otherwise
            
            y_relative = 0.013583333333333;
            abs_length = y_relative*((ylim_custom(2)-ylim_custom(1))/ax_position(4));

            % 在原有刻度上添加额外的标记
            for i = 1:length(positions)
                x_position = positions(i);
                x_norm2 = (x_position - xlim_custom(1)) / (xlim_custom(2) - xlim_custom(1));

                x_relative = x_norm2 * ax_position(3);

                % 绘制垂直于 x 轴的线条
                line([x_relative*(26/ax_position(3)), x_relative*(26/ax_position(3))], [abs_length*1.3 + ylim_custom(1),ylim_custom(1)], ...
                    'Color', 'k', 'LineStyle', '-', 'LineWidth', 3, 'HandleVisibility', 'off');
            end
            
            % 调整 x 轴标签位置
            y_relative3 = 0.040750000000000;
            abs_length3 = y_relative3*((ylim_custom(2)-ylim_custom(1))/ax_position(4));

            ax.XAxis.Label.Position(2) = ax.XAxis.Label.Position(2) - abs_length3;          
    end

    % 循环添加标签
    for i = 1:length(positions)
        x_label = positions(i);
        x_norm_label = (x_label - xlim_custom(1)) / (xlim_custom(2) - xlim_custom(1));

        x_relative_label = x_norm_label * ax_position(3);

        abs_length = y_relative*((ylim_custom(2)-ylim_custom(1))/ax_position(4));
        
        % 添加标签
        text(x_relative_label*(26/ax_position(3)), ylim_custom(1)-abs_length, labels{i}, 'Color', 'k', 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 20,'FontWeight', 'bold');
    end
    
    box off     % 取消边框
    ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
    set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度
    
    % 设置 x 轴和 y 轴的线条粗细
    ax1.XAxis.LineWidth = 2;
    ax1.YAxis.LineWidth = 2;
    ax1.XAxis.FontSize = 20;
    ax1.YAxis.FontSize = 20;
    ax1.YAxis.FontWeight = 'bold';
    
    % 设置保存路径和文件名
    save_path = 'E:\Data\站点数据_1217\绘图数据\Pictures';

    % 确保保存路径存在，如果不存在则创建
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    
    % 完整的保存路径和文件名
    full_file_path = fullfile(save_path, output_file);
    
    % 保存图形为 JPEG 格式，设置分辨率
    saveas(gcf, full_file_path, 'svg');
end
