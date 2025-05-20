clc;
clear;

% 自定义 x 轴标签
% 定义要标注的位置和对应的标签
positions = [1, 5, 9, 12, 14, 17, 21, 25];
labels = {'-3M', '-2M', '-1M', '-8D', '8D', '1M', '2M', '3M'};

% 设置默认字体为 Times New Roman
set(groot, 'defaultAxesFontName', 'Times New Roman');

% 指定路径
path = 'E:\Data\Drought_response\Season\FD\1CI';

% 定义四个季节的名称
seasons = {'spring', 'summer', 'autumn', 'winter'};

% 初始化一个 cell 数组来存储季节文件夹的完整路径
season_paths = cell(size(seasons));

% 构建每个季节文件夹的完整路径
for i = 1:numel(seasons)
    season_paths{i} = fullfile(path, seasons{i});
end

% 现在您可以在每个季节文件夹路径上运行之前修改过的处理代码
for season_idx = 1:numel(season_paths)
    season_folder = season_paths{season_idx};

    [~, season_file] = fileparts(season_folder);

    season_file(1) = upper(season_file(1));  % 将第一个字符转换为大写

    % 使用 filesep 获取文件分隔符（根据操作系统不同可能是 '\' 或 '/'）
    filesep_char = filesep;
    
    % 使用 split 函数拆分路径字符串
    folders = split(season_folder, filesep_char);
    
    % 获取指定路径下的所有文件和文件夹
    contents = dir(season_folder);
    
    % 初始化一个额外处理的 cell 数组来存储文件名
    selected_files_periods = {};
    
    selected_files_count = {};

    % 遍历所有内容
    for i = 1:numel(contents)
        % 忽略当前目录（'.'）和父目录（'..'）
        if strcmp(contents(i).name, '.') || strcmp(contents(i).name, '..')
            continue;
        end
        
        % 如果当前内容是文件夹，则进一步处理
        if contents(i).isdir
            folder_name = fullfile(season_folder, contents(i).name);
            disp(['读取处理 文 件 夹 ：', folder_name]);
    
            % 获取子文件夹名称
            [~, folder_name_short, ~] = fileparts(folder_name);
            
            % 获取当前文件夹中的所有文件和文件夹
            sub_contents = dir(folder_name);
    
            % 初始化一个 cell 数组来存储文件名
            selected_files = {};
            
            % 遍历当前文件夹中的所有内容
            for j = 1:numel(sub_contents)
                % 忽略当前目录（'.'）和父目录（'..'）
                if strcmp(sub_contents(j).name, '.') || strcmp(sub_contents(j).name, '..')
                    continue;
                end
                        
                % 如果当前内容是文件且文件名包含 "mean" 字符串且不包含 "period" 字符串
                if ~sub_contents(j).isdir && contains(sub_contents(j).name, 'mean') && ~contains(sub_contents(j).name, 'period')
                    file_name = fullfile(folder_name, sub_contents(j).name);
                    
                    % 将符合条件的文件路径添加到 selected_files 数组中
                    selected_files{end+1} = file_name;
                end
                % 如果 folder_name_short 是 'Mean'，进行额外处理
                if strcmp(folder_name_short, 'Mean')
    
                    if ~sub_contents(j).isdir && contains(sub_contents(j).name, 'mean') && contains(sub_contents(j).name, 'period')
                        file_name = fullfile(folder_name, sub_contents(j).name);
                        selected_files_periods{end+1} = file_name;
                    end
                    if ~sub_contents(j).isdir && contains(sub_contents(j).name, 'count')
                        file_name = fullfile(folder_name, sub_contents(j).name);
                        selected_files_count{end+1} = file_name;
                    end
                end
            end
    
            % 动态生成变量名并存储文件列表
            eval(['selected_files_', folder_name_short, '= selected_files;']);
        end
    end
    
    for m = 1:numel(selected_files_Fast)
        file_name_fast = selected_files_Fast{m};
        dataArray_fast = readmatrix(file_name_fast);
        dataArray_fast = dataArray_fast(:)'; 
    
        file_name_mean = selected_files_Mean{m};
        dataArray_mean = readmatrix(file_name_mean);
        dataArray_mean = dataArray_mean(:)'; 
    
        file_name_slow = selected_files_Slow{m};
        dataArray_slow = readmatrix(file_name_slow);
        dataArray_slow = dataArray_slow(:)';
    
        % 提取文件名部分（不包括路径和扩展名）
        [~, base_name, ~] = fileparts(file_name_fast);
        words = split(base_name, '_');
        % 获取第三个和第四个单词
        third_word = words{3};
        fourth_word = words{4};
    
        % 构建新的输出文件名
        output_file = sprintf('%s_%s_%s.svg', fourth_word, third_word, season_file);
    
        % 根据 fourth_word 选择颜色
        switch fourth_word
            case 'ESI'
%                 color_fast = '#F12C32';
                color_mean = '#F12C32';
%                 color_slow = '#FBA250'; 
            case 'GPP'
                color_fast = '#112F97';
                color_mean = '#095BD0';
                color_slow = '#87BFDF';
            case 'NPP'
                color_fast = '#691E96';
                color_mean = '#9A68C4';
                color_slow = '#D6BFE2';
            case 'T'
                color_fast = '#0C721F';
                color_mean = '#61AF25';
                color_slow = '#B5DC6C';
            case 'NEP'
                color_fast = '#11AB0B';
                color_mean = '#17D80F';
                color_slow = '#72F46D';
            case 'Ra'
                color_fast = '#FAA103';
                color_mean = '#FFD600';
                color_slow = '#FFF176';
            case 'Reco'
                color_fast = '#D25203';
                color_mean = '#FF780D';
                color_slow = '#FCA16B';
            case 'LAI'
                color_fast = '#27A09E';
                color_mean = '#30CE88';
                color_slow = '#7DE393';  
            case 'SIF'
                color_fast = '#4FFD95';
                color_mean = '#97FD52';
                color_slow = '#C5FB2E';
        end
    
        period = 1:25;
        serious_period = 13;
    
        % 创建新的图窗并设置其位置和大小
        figure('Position', [800, 100, 800, 580]);
        
        % 绘制折线图
%         plot(period, dataArray_fast, '-','Color', color_fast, 'LineWidth', 3);
        hold on;
        plot(period, dataArray_mean, '-','Color', color_mean, 'LineWidth', 4);
%         plot(period, dataArray_slow, '-','Color', color_slow, 'LineWidth', 3);
    
        period_filename = selected_files_periods{m};
        period_points = readmatrix(period_filename);
        period_points = period_points(:)';
    
        count_filename = selected_files_count{m};
        drought_count = readmatrix(count_filename);
        drought_number = drought_count(13);
    
        % 根据 fourth_word 设置不同的 x 和 y 轴范围和刻度间隔
        switch fourth_word
            case 'ESI'
                xlim_custom = [0 26];        
                xticks_custom = 0:1:26;      
                ylim_custom = [-1.5 1.5];    
                yticks_custom = -1.5:0.5:1.5;
            case 'GPP'
                xlim_custom = [0 26];        
                xticks_custom = 0:1:26;      
                ylim_custom = [-2 14];    
                yticks_custom = -2:2:14;
            case 'NPP'
                xlim_custom = [0 26];
                xticks_custom = 0:1:26;
                ylim_custom = [-1 6];
                yticks_custom = -1:1:6;
            case 'T'
                xlim_custom = [0 26];        
                xticks_custom = 0:1:26;      
                ylim_custom = [-0.8 4.8];    
                yticks_custom = -0.8:0.8:4.8;
            case 'NEP'
                xlim_custom = [0 26];
                xticks_custom = 0:1:26;
                ylim_custom = [-1.5 2];
                yticks_custom = -1.5:0.5:2;
            case 'Ra'
                xlim_custom = [0 26];
                xticks_custom = 0:1:26;
                ylim_custom = [-2 10];
                yticks_custom = -2:2:10;
            case 'Reco'
                xlim_custom = [0 26];
                xticks_custom = 0:1:26;
                ylim_custom = [-2 12];
                yticks_custom = -2:2:12;
            case 'LAI'
                xlim_custom = [0 26];        
                xticks_custom = 0:1:26;      
                ylim_custom = [0 6];    
                yticks_custom = 0:1:6;
            case 'SIF'
                xlim_custom = [0 26];        
                xticks_custom = 0:1:26;      
                ylim_custom = [-0.2 1];    
                yticks_custom = -0.2:0.2:1;
        end
    
        % 设置 x 和 y 轴范围
        xlim(xlim_custom);
        ylim(ylim_custom);
        
        % 设置 x 和 y 轴刻度间隔
        xticks(xticks_custom);
        yticks(yticks_custom);
    
        switch fourth_word
            case 'ESI'
                % 设置 y 轴刻度，只显示部分刻度
                set(gca, 'YTick', [-1.0, -0.5, 0, 0.5, 1]);
            case 'GPP'
                set(gca, 'YTick', [0, 3, 6, 9, 12]);
            case 'NPP'
                set(gca, 'YTick', [0, 1, 2, 3, 4, 5]);
            case 'T'
                set(gca, 'YTick', [0, 0.8, 1.6, 2.4, 3.2, 4.0]);
            case 'NEP'
                set(gca, 'YTick', [-1.0, -0.5, 0, 0.5, 1.0, 1.5]);
            case 'Ra'
                set(gca, 'YTick', [0, 2, 4, 6, 8]);
            case 'Reco'
                set(gca, 'YTick', [0, 2, 4, 6, 8, 10]);
            case 'LAI'
                set(gca, 'YTick', [1, 2, 3, 4, 5]);
            case 'SIF'
                set(gca, 'YTick', [0, 0.2, 0.4, 0.6, 0.8]);
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
%         plot(period, dataArray_fast, '-','Color', color_fast, 'LineWidth', 3);
        plot(period, dataArray_mean, '-','Color', color_mean, 'LineWidth', 4);
%         plot(period, dataArray_slow, '-','Color', color_slow, 'LineWidth', 3);
    
        % 设置标题和 y 轴标签，根据 third_word 和 drought_speed
        if strcmp(third_word, 'FD')
            if strcmp(fourth_word, 'ESI')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = '';
            elseif strcmp(fourth_word, 'GPP')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = 'g C m⁻² d⁻¹';
            elseif strcmp(fourth_word, 'NPP')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = 'g C m⁻² d⁻¹';
            elseif strcmp(fourth_word, 'T')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = 'mm day⁻¹';
            elseif strcmp(fourth_word, 'NEP')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = 'g C m⁻² d⁻¹';
            elseif strcmp(fourth_word, 'Ra')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = 'g C m⁻² d⁻¹';
            elseif strcmp(fourth_word, 'Reco')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = 'g C m⁻² d⁻¹';
            elseif strcmp(fourth_word, 'LAI')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = '';
            elseif strcmp(fourth_word, 'SIF')
                title_str = sprintf('%s Response to %s Flash Drought', fourth_word, season_file);
                ylabel_str = 'mWm⁻² nm⁻¹ sr⁻¹';
            end
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
        ylabel(ylabel_str, 'FontName', 'Times New Roman', 'FontSize', 24,'FontWeight', 'bold');
        xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');

        if strcmp(third_word, 'FD') 
            legend_handle = legend('Average', 'Location', 'northeast','FontSize', 22, 'FontWeight', 'bold');
        elseif strcmp(third_word, 'SD') 
            legend_handle = legend('Fastest','Average','Slowest', 'Location', 'northeast','FontSize', 22, 'FontWeight', 'bold');
        end
        legend_handle.Box = 'off';

        ax.YAxis.Label.Position(1) = -2;
    
        % 根据不同的 fourth_word 来调整添加额外的垂直线条和调整 x 轴标签位置
        switch fourth_word
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

        if strcmp(third_word, 'FD') 
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
            y_distance = ax_position(2) + ax_position(4) -distance_top_boundary_min-legend_position(4);

            % 在相对坐标系中添加文本
            annotation_handle = annotation('textbox', [x_distance, y_distance, 0.1, 0.1], ...
                'String', sprintf('N:%d', drought_number), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                'FontSize', 20, 'FontWeight', 'bold','FontName', 'Times New Roman', 'EdgeColor', 'none','FitBoxToText', 'on');
        elseif strcmp(third_word, 'SD') 
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
            y_distance = ax_position(2) + ax_position(4) -distance_top_boundary_min-legend_position(4)/3;

            % 在相对坐标系中添加文本
            annotation_handle = annotation('textbox', [x_distance, y_distance, 0.1, 0.1], ...
                'String', sprintf('N:%d', drought_number), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'EdgeColor', 'none','FitBoxToText', 'on');
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
        save_path = 'E:\Data\Drought_quantification\response_0417';
    
        % 确保保存路径存在，如果不存在则创建
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        
        % 完整的保存路径和文件名
        full_file_path = fullfile(save_path, output_file);
        
        % 设置输出分辨率为 500 DPI
        resolution = 800;
        
        % 保存图形为 JPEG 格式，设置分辨率
        saveas(gcf, full_file_path, 'svg');
    
    end
end
