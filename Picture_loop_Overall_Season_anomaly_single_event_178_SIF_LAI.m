clc;
clear;

% 自定义 x 轴标签
% 定义要标注的位置和对应的标签
positions = [1, 5, 9, 12, 14, 17, 21, 25];
labels = {'-3M', '-2M', '-1M', '-8D', '8D', '1M', '2M', '3M'};

% 设置默认字体为 Times New Roman
set(groot, 'defaultAxesFontName', 'Times New Roman');

% 指定路径
path = 'E:\Data\站点数据_1217\千烟洲\FD_178_2';

% 获取当前文件夹中的所有文件和文件夹
sub_contents = dir(path);
 
% 初始化一个 cell 数组来存储文件名
selected_files = {};
selected_files_periods = {};
selected_files_count = {};

% 遍历当前文件夹中的所有内容
for j = 1:numel(sub_contents)
    % 忽略当前目录（'.'）和父目录（'..'）
    if strcmp(sub_contents(j).name, '.') || strcmp(sub_contents(j).name, '..')
        continue;
    end
            
    % 如果当前内容是文件且文件名包含 "mean" 字符串且不包含 "period" 字符串
    if ~sub_contents(j).isdir && contains(sub_contents(j).name, 'mean') && ~contains(sub_contents(j).name, 'period')
        file_name = fullfile(path, sub_contents(j).name);
        
        % 将符合条件的文件路径添加到 selected_files 数组中
        selected_files{end+1} = file_name;
    end

    if ~sub_contents(j).isdir && contains(sub_contents(j).name, 'period')
        file_name = fullfile(path, sub_contents(j).name);
        selected_files_periods{end+1} = file_name;
    end
    if ~sub_contents(j).isdir && contains(sub_contents(j).name, 'count')
        file_name = fullfile(path, sub_contents(j).name);
        selected_files_count{end+1} = file_name;
    end
end

% % 动态生成变量名并存储文件列表
% eval(['selected_files_', folder_name_short, '= selected_files;']);
% eval(['selected_files_periods_', folder_name_short, '= selected_files_periods;']);
% eval(['selected_files_count_', folder_name_short, '= selected_files_count;']);
% 
% 将这些变量放入一个 cell 数组中
selected_files_all = {selected_files};

for i = 1:length(selected_files_all)
    current_files = selected_files_all{i};

    % 创建新的图窗并设置其位置和大小
    figure('Position', [1700, 100, 800, 580]);
    hold on;
    for j = 1:length(current_files)
        % 提取当前文件
        current_file = current_files{j};

        [~, site_name] = fileparts(fileparts(fileparts(current_file)));

        % 提取文件名部分（不包括路径和扩展名）
        [~, base_name, ~] = fileparts(current_file);
        words = split(base_name, '_');

        % 获取第三个和第四个单词
        third_word = words{3};
        fourth_word = words{4};

        % 构建新的输出文件名
        output_file = sprintf('%s_%s_%s_only_4.svg', words{1}, words{2}, site_name);

        period = 1:25;
        serious_period = 13;

        % 根据 fourth_word
        switch fourth_word
%             case 'ESI'
%                 data_ESI = readmatrix(current_file);
%                 data_ESI = data_ESI(:)';
% 
%                 % 对数据进行距平处理
%                 mean_value_ESI = mean(data_ESI);
%                 anomaly_ESI = data_ESI - mean_value_ESI;
% 
%                 plot(period, anomaly_ESI, '-','Color', '#ED2F33', 'LineWidth', 3);
            case 'LAI'
                data_LAI = readmatrix(current_file);
                data_LAI = data_LAI(:)';

                yyaxis left % 激活左侧 y 轴
                plot(period, data_LAI, '-','Color', '#85EDB7', 'LineWidth', 3);
            case 'SIF'
                data_SIF = readmatrix(current_file);
                data_SIF = data_SIF(:)';

                yyaxis right % 激活右侧 y 轴
                plot(period, data_SIF, '-','Color', '#C5FB2E', 'LineWidth', 3);
%             case 'T'
%                 data_T = readmatrix(current_file);
%                 data_T = data_T(:)';
% 
%                 yyaxis right % 激活右侧 y 轴
%                 plot(period, data_T, '-','Color', '#0C721F', 'LineWidth', 3);
%             case 'NEP'
%                 data_NEP = readmatrix(current_file);
%                 data_NEP = data_NEP(:)';
% 
%                 yyaxis left % 激活左侧 y 轴
%                 plot(period, data_NEP, '-','Color', '#17D80F', 'LineWidth', 3);
%             case 'GPP'
%                 data_GPP = readmatrix(current_file);
%                 data_GPP = data_GPP(:)';
% 
%                 yyaxis left % 激活左侧 y 轴
%                 plot(period, data_GPP, '-', 'Color', [9/255, 91/255, 208/255], 'LineWidth', 3);
        end
    end

    periods_files = selected_files_periods;

    % 读取第一个文件的数据作为基准
    first_filename = periods_files{1};
    period_points = readmatrix(first_filename);
    period_points = period_points(:)'; % 转换为行向量

    count_filename = selected_files_count{1};
    drought_count = readmatrix(count_filename);
    drought_number = drought_count(13);

    % 定义标签字符串
%     ylabel_left_str = '( GPP/NEP ) g C m⁻² d⁻¹';
%     ylabel_right_str = '( T ) mm day⁻¹';

    ylabel_left_str = 'LAI';
    ylabel_right_str = 'SIF (mWm⁻² nm⁻¹ sr⁻¹)';

    yyaxis right

    h_right_label = ylabel(ylabel_right_str, 'FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold');

    xlim_custom = [0 26];
    xticks_custom = 0:1:26;
    ylim_custom = [-0.05 0.65];
    yticks_custom = -0.05:0.1:0.65;

    % 设置 x 和 y 轴范围
    xlim(xlim_custom);
    ylim(ylim_custom);
    
    % 设置 x 和 y 轴刻度间隔
    xticks(xticks_custom);
    yticks(yticks_custom);
    
    % 设置 y 轴刻度，只显示部分刻度
    set(gca, 'YTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]);
    
    % 获取当前的 y 轴刻度
    yTicks = get(gca, 'YTick');
    
    % 将刻度四舍五入到一位小数，并转换成字符串格式
    formattedLabels = sprintf('%.1f\n', yTicks);
    
    % 设置 y 轴刻度的文本格式
    set(gca, 'YTickLabel', formattedLabels);
    set(gca, 'YColor', [0, 0, 0]);

    yyaxis left

    h_left_label = ylabel(ylabel_left_str, 'FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold');

    xlim_custom = [0 26];        
    xticks_custom = 0:1:26;      
    ylim_custom = [-0.5 5.5];    
    yticks_custom = -0.5:1:5.5;

    % 设置 x 和 y 轴范围
    xlim(xlim_custom);
    ylim(ylim_custom);
    
    % 设置 x 和 y 轴刻度间隔
    xticks(xticks_custom);
    yticks(yticks_custom);
    
    % 设置 y 轴刻度，只显示部分刻度
    set(gca, 'YTick', [0, 1, 2, 3, 4, 5]);
    
    % 获取当前的 y 轴刻度
    yTicks = get(gca, 'YTick');
    
    % 将刻度四舍五入到一位小数，并转换成字符串格式
    formattedLabels = sprintf('%.1f\n', yTicks);
    
    % 设置 y 轴刻度的文本格式
    set(gca, 'YTickLabel', formattedLabels);
    set(gca, 'YColor', [0, 0, 0]);

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

    for j = 1:length(current_files)
        % 提取当前文件
        current_file = current_files{j};

        [~, site_name] = fileparts(fileparts(current_file));

        % 提取文件名部分（不包括路径和扩展名）
        [~, base_name, ~] = fileparts(current_file);
        words = split(base_name, '_');

        % 获取第三个和第四个单词
        third_word = words{3};
        fourth_word = words{4};

        % 根据 fourth_word
        switch fourth_word
%             case 'ESI'
%                 data_ESI = readmatrix(current_file);
%                 data_ESI = data_ESI(:)';
% 
%                 % 对数据进行距平处理
%                 mean_value_ESI = mean(data_ESI);
%                 anomaly_ESI = data_ESI - mean_value_ESI;
% 
%                 h_ESI = plot(period, anomaly_ESI, '-','Color', '#ED2F33', 'LineWidth', 3);
            case 'LAI'
                data_LAI = readmatrix(current_file);
                data_LAI = data_LAI(:)';

                yyaxis left % 激活左侧 y 轴
                h_LAI = plot(period, data_LAI, '-','Color', '#85EDB7', 'LineWidth', 3);

                y_min = min(ylim); % 获取 y 轴范围的最小值
                y_max = max(ylim); % 获取 y 轴范围的最大值
                line([0.12, 0.12], [y_min, y_max], 'Color', '#85EDB7', 'LineWidth', 3);
            case 'sif'
                data_SIF = readmatrix(current_file);
                data_SIF = data_SIF(:)';

                yyaxis right % 激活左侧 y 轴
                h_SIF = plot(period, data_SIF, '-','Color', '#C5FB2E', 'LineWidth', 3);

                y_min = min(ylim); % 获取 y 轴范围的最小值
                y_max = max(ylim); % 获取 y 轴范围的最大值
                line([25.88, 25.88], [y_min, y_max], 'Color', '#C5FB2E', 'LineWidth', 3);
%             case 'T'
%                 data_T = readmatrix(current_file);
%                 data_T = data_T(:)';
% 
%                 yyaxis right % 激活右侧 y 轴
%                 h_T = plot(period, data_T, '-','Color', '#0C721F', 'LineWidth', 3);
% 
%                 y_min = min(ylim); % 获取 y 轴范围的最小值
%                 y_max = max(ylim); % 获取 y 轴范围的最大值
%                 line([25.86, 25.86], [y_min, y_max], 'Color', '#0C721F', 'LineWidth', 4);
%             case 'NEP'
%                 data_NEP = readmatrix(current_file);
%                 data_NEP = data_NEP(:)';
% 
%                 yyaxis left % 激活左侧 y 轴
%                 h_NEP = plot(period, data_NEP, '-','Color', '#17D80F', 'LineWidth', 3);
%             case 'GPP'
%                 data_GPP = readmatrix(current_file);
%                 data_GPP = data_GPP(:)';
% 
%                 yyaxis left % 激活左侧 y 轴
%                 h_GPP = plot(period, data_GPP, '-', 'Color', [9/255, 91/255, 208/255], 'LineWidth', 3);
        end
    end

    % 设置标题和 y 轴标签，根据 third_word 和 drought_speed
    title_str = sprintf('Response to a Single Flash Drought');
    ylabel_str = '';

    % 获取当前坐标轴对象
    ax = gca;
    
    % 设置 x 轴线条粗细
    ax.XAxis.LineWidth = 2;
    ax.XAxis.FontSize = 20;

    yyaxis right
    ax = gca;
    ax.YAxis(2).TickLength = [0.015; 0];
    ax.YAxis(2).LineWidth = 2; % 左侧 y 轴的线宽
    ax.YAxis(2).FontSize = 20; % 左侧 y 轴的字体大小
    ax.YAxis(2).FontWeight = 'bold'; % 左侧 y 轴字体加粗

    h_right_label.Position(1) = 27.5;

    set(gca,  'XTickLabel', []);

    % 激活左侧 y 轴并设置属性
    yyaxis left
    ax = gca;
    ax.YAxis(1).TickLength = [0.015; 0];
    ax.YAxis(1).LineWidth = 2; % 左侧 y 轴的线宽
    ax.YAxis(1).FontSize = 20; % 左侧 y 轴的字体大小
    ax.YAxis(1).FontWeight = 'bold'; % 左侧 y 轴字体加粗

    h_left_label.Position(1) = -2; % 左右方向调整

%     set(gca,  'XTickLabel', []);

    % 获取坐标轴位置信息
    ax_position = ax.Position;

    % 设置标题
    title_handle_1 = title(title_str, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');

    % 设置标签
    xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold');

%     legend_handle = legend([h_GPP, h_NEP, h_T], ...
%     {'GPP', 'NEP', 'T'}, ...
%     'Location', 'northeast', 'FontSize', 22, 'FontWeight', 'bold');

    legend_handle = legend([h_SIF, h_LAI], ...
    {'SIF', 'LAI'}, ...
    'Location', 'northeast', 'FontSize', 22, 'FontWeight', 'bold');

    legend_handle.Box = 'off';

    y_norm2 = 0.1 / (ylim_custom(2) - ylim_custom(1));
%     y_relative = y_norm2 * ax_position(4);
    y_relative = 0.013583333333333;
    abs_length = y_relative*((ylim_custom(2)-ylim_custom(1))/ax_position(4));

    % 在原有刻度上添加额外的标记
    for y = 1:length(positions)
        x_position = positions(y);
        x_norm2 = (x_position - xlim_custom(1)) / (xlim_custom(2) - xlim_custom(1));

        x_relative = x_norm2 * ax_position(3);

        % 绘制垂直于 x 轴的线条
        line([x_relative*(26/ax_position(3)), x_relative*(26/ax_position(3))], [abs_length*1.3 + ylim_custom(1),ylim_custom(1)], ...
            'Color', 'k', 'LineStyle', '-', 'LineWidth', 3, 'HandleVisibility', 'off');
    end
    
    % 调整 x 轴标签位置
    y_norm3 = 0.3 / (ylim_custom(2) - ylim_custom(1));
%     y_relative3 = y_norm3 * ax_position(4);
    y_relative3 = 0.040750000000000;
    abs_length3 = y_relative3*((ylim_custom(2)-ylim_custom(1))/ax_position(4));

    ax.XAxis.Label.Position(2) = ax.XAxis.Label.Position(2) - abs_length3;

    % 循环添加标签
    for o = 1:length(positions)
        x_label = positions(o);
        x_norm_label = (x_label - xlim_custom(1)) / (xlim_custom(2) - xlim_custom(1));

        x_relative_label = x_norm_label * ax_position(3);

        abs_length = y_relative*((ylim_custom(2)-ylim_custom(1))/ax_position(4));
        
        % 添加标签
        text(x_relative_label*(26/ax_position(3)), ylim_custom(1)-abs_length, labels{o}, 'Color', 'k', 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 20,'FontWeight', 'bold');
    end

    % 获取图例的位置信息
    legend_position = get(legend_handle, 'Position');

    % 计算图例的右边界和上边界位置
    legend_right = legend_position(1) + legend_position(3);
    legend_top = legend_position(2) + legend_position(4);
    
    % 计算图例与 xlim_custom 和 ylim_custom 边界的距离
    distance_right_boundary_min = ax_position(1) + ax_position(3) - legend_right;
    distance_top_boundary_min = ax_position(2) + ax_position(4) - legend_top;
  
    % 计算文本位置
    x_distance = ax_position(1) + distance_right_boundary_min;
    y_distance = ax_position(2) + ax_position(4) -distance_top_boundary_min-legend_position(4)/2;

    % 在相对坐标系中添加文本
    annotation_handle = annotation('textbox', [x_distance, y_distance, 0.1, 0.1], ...
        'String', sprintf('N:%d', drought_number), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'FontSize', 20, 'FontWeight', 'bold','FontName', 'Times New Roman', 'EdgeColor', 'none','FitBoxToText', 'on');

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
