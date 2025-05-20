clc;
clear;

% 自定义 x 轴标签
% 定义要标注的位置和对应的标签
positions = [1, 5, 9, 12, 14, 17, 21, 25];
labels = {'-3M', '-2M', '-1M', '-8D', '8D', '1M', '2M', '3M'};

% 设置默认字体为 Times New Roman
set(groot, 'defaultAxesFontName', 'Times New Roman');

period = 1:25;

serious_period = 13;

% 指定路径
path = 'E:\Data\Drought_quantification\Drought_quantify\Season_ALL_1_month\FD\1CI_new\1CI';
% path = 'E:\Data\Drought_quantification\Drought_quantify\Season_N&P_1_month\Natural\FD\1CI_new\1CI';
% path = 'E:\Data\Drought_quantification\Drought_quantify\Season_N&P_1_month\Planted\FD\1CI_new\1CI';

[~,segmented_name] = fileparts(fileparts(fileparts(fileparts(path))));

if contains(segmented_name, 'ALL')
    segmented_name = 'ALL';
end

% 定义四个季节的名称
seasons = {'spring', 'summer', 'autumn', 'winter'};

% 初始化一个 cell 数组来存储季节文件夹的完整路径
season_paths = cell(size(seasons));

% 构建每个季节文件夹的完整路径
for i = 1:numel(seasons)
    season_paths{i} = fullfile(path, seasons{i});
end

% 现在您可以在每个季节文件夹路径上运行之前修改过的处理代码
for season_idx = 3:numel(season_paths) 
    season_folder = season_paths{season_idx};

    [~, season_file] = fileparts(season_folder);

    season_file(1) = upper(season_file(1));  % 将第一个字符转换为大写

    % 使用 filesep 获取文件分隔符（根据操作系统不同可能是 '\' 或 '/'）
    filesep_char = filesep;
    
    % 使用 split 函数拆分路径字符串
    folders = split(season_folder, filesep_char);
    
    % 获取指定路径下的所有文件和文件夹
    contents = dir(season_folder); 
    
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

            % 检查子文件夹是否为 'Mean'
            if strcmp(folder_name_short, 'Mean')
            
                % 获取当前文件夹中的所有文件和文件夹
                sub_contents = dir(folder_name);

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
                        file_name = fullfile(folder_name, sub_contents(j).name);

                        % 提取文件名部分（不包括路径和扩展名）
                        [~, base_name, ~] = fileparts(file_name);
                        words = split(base_name, '_');
                        % 获取第四个单词
                        fourth_word = words{4};

                        if ismember(fourth_word, {'T','NEP', 'GPP', 'NPP', 'Ra', 'Reco'})
                            % 将符合条件的文件路径添加到 selected_files 数组中
                            selected_files{end+1} = file_name;
                        end
                    end
        
                    if ~sub_contents(j).isdir && contains(sub_contents(j).name, 'mean') && contains(sub_contents(j).name, 'period')
                        file_name = fullfile(folder_name, sub_contents(j).name);

                        % 提取文件名部分（不包括路径和扩展名）
                        [~, base_name, ~] = fileparts(file_name);
                        words = split(base_name, '_');
                        % 获取第四个单词
                        fourth_word = words{4};

                        if ismember(fourth_word, {'T','NEP', 'GPP', 'NPP', 'Ra', 'Reco'})
                            selected_files_periods{end+1} = file_name;
                        end
                    end
                    if ~sub_contents(j).isdir && contains(sub_contents(j).name, 'count')
                        file_name = fullfile(folder_name, sub_contents(j).name);

                        % 提取文件名部分（不包括路径和扩展名）
                        [~, base_name, ~] = fileparts(file_name);
                        words = split(base_name, '_');
                        % 获取第四个单词
                        fourth_word = words{4};
                        
                        if ismember(fourth_word, {'T','NEP', 'GPP', 'NPP', 'Ra', 'Reco'})
                            selected_files_count{end+1} = file_name;
                        end
                    end
                end
            end
        end
    end

    for m = 1:numel(selected_files)
        selected_file = selected_files{m};
        [~, tree_kind] = fileparts(fileparts(fileparts(fileparts(fileparts(fileparts(fileparts(selected_file)))))));

        if strcmp(tree_kind, 'Natural') || strcmp(tree_kind, 'Planted')
            % tree_kind 保持不变
        elseif contains(tree_kind, 'ALL')
            tree_kind = '';
        end

        % 提取文件名部分（不包括路径和扩展名）
        [~, data_kinds, ~] = fileparts(selected_file);
        data_words = split(data_kinds, '_');
        % 获取第四个单词
        third_word = data_words{3};
        data_kind = data_words{4};

        if strcmp(data_kind, 'GPP')
            data_GPP = readmatrix(selected_file);
            data_GPP = data_GPP(:)';

            % 获取 period_points
            for n = 1:numel(selected_files_periods)
                selected_file_period = selected_files_periods{n};
                [~, period_kinds, ~] = fileparts(selected_file_period);
                period_words = split(period_kinds, '_');
                period_kind = period_words{4};
    
                if strcmp(period_kind, 'GPP') % 假设period_kind也应该是GPP
                    period_points = readmatrix(selected_file_period);
                    period_points = period_points(:)';
                    break; % 找到匹配的文件就退出循环
                end
            end

            Interpolation_values_drought = interp1(period, data_GPP, period_points);

            selected_file_path_parts = split(selected_file, filesep);
            fourth_file_name = selected_file_path_parts(5);

            if strcmp(fourth_file_name, 'Season_ALL_1_month')
                % 构造对应文件的路径
                corresponding_file = strrep(selected_file, 'Season_ALL_1_month', 'Season_ALL_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(10) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s.svg', data_kind, third_word, season_file);
            else
                corresponding_file = strrep(selected_file, 'Season_N&P_1_month', 'Season_N&P_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(11) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s(%s).svg', data_kind, third_word, season_file, tree_kind);
            end
            
            % 重新组合路径
            corresponding_file = fullfile(path_parts{:});

            data_GPP_no_drought = readmatrix(corresponding_file);
            data_GPP_no_drought = data_GPP_no_drought(:)';

            Interpolation_values_no_drought = interp1(period, data_GPP_no_drought, period_points);

            % 提取从 start_period 开始的时间和蒸腾变化量数据
            idx_start = find(period >= ceil(period_points(1)), 1);
            idx_end = find(period >= ceil(period_points(2)), 1);
            drought_period = [period_points(1), period(idx_start:idx_end-1),period_points(2)];
            drought_GPP_rate = [Interpolation_values_drought(1), data_GPP(idx_start:idx_end-1),Interpolation_values_drought(2)];
            
            % 干旱期间的累计蒸腾量
            cumulative_transpiration_droughting = cumtrapz(drought_period, drought_GPP_rate);
            
            new_drought_period = [drought_period, period(idx_end:end)];
            new_drought_GPP_rate = [drought_GPP_rate, data_GPP(idx_end:end)];
            
            cumulative_transpiration_drought = cumtrapz(new_drought_period, new_drought_GPP_rate);
            
            periods_before_drought = [period(1:idx_start - 1), period_points(1)];
            GPP_rate_before_drought = [data_GPP(1:idx_start - 1), Interpolation_values_drought(1)];
            cumulative_transpiration_before_drought = cumtrapz(periods_before_drought, GPP_rate_before_drought);
            
            new_period = [period(1:idx_start - 1), drought_period, period(idx_end:end)];
            new_GPP_rate = [data_GPP(1:idx_start - 1), drought_GPP_rate, data_GPP(idx_end:end)];
            
            % 假设没有发生干旱 假设没有发生干旱 假设没有发生干旱
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况

            integer_periods = ceil(period_points(1)):floor(period_points(2));
            
            % 计算整个时间范围内的累积蒸腾变化量
            cumulative_transpiration_total = cumtrapz(new_period, new_GPP_rate);
            
            lush_period = [period_points(1), integer_periods, period_points(2)];
            lush_GPP_rate = [Interpolation_values_no_drought(1), data_GPP_no_drought(integer_periods), Interpolation_values_no_drought(2)]; 
            
            % 计算没有干旱发生 的累计蒸腾量 针对干旱时期对应的时间
            cumulative_transpiration_lush = cumtrapz(lush_period, lush_GPP_rate);
            
            new_lush_period = [lush_period, period(idx_end:end)];
            new_lush_GPP_rate = [lush_GPP_rate, data_GPP_no_drought(idx_end:end)];
           
            % 计算没有干旱发生 的累计蒸腾量 针对干旱开始到25！！
            cumulative_transpiration_lush_after_drought = cumtrapz(new_lush_period, new_lush_GPP_rate);
            
            new_lush_period_total = [period(1:idx_start - 1), new_lush_period];
            new_lush_GPP_rate_total = [data_GPP_no_drought(1:idx_start - 1), new_lush_GPP_rate];
            
            cumulative_transpiration_lush_total = cumtrapz(new_lush_period_total, new_lush_GPP_rate_total);
            
            % 创建图形窗口并调整尺寸
            % figure('Units', 'normalized', 'Position', [0.6, 0.1, 0.3, 0.5]);
            figure('Units', 'pixels', 'Position', [750,100,800,1000]);
            % 绘制蒸腾变化量随时间变化的曲线及其填充
            h1 = subplot(2, 1, 1);
            % 绘制干旱的情况
            p_star = plot(drought_period([1, end]), drought_GPP_rate([1, end]), 'p', 'MarkerSize', 3, 'LineWidth', 2 ,...
                'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255]);
            hold on;
%             p_original = plot(period, data_GPP, '-','Color', [9/255, 91/255, 208/255, 0.8],'LineWidth', 1.5);
%   
%             p_drought = plot(drought_period, drought_GPP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.5); % 底部边界线
%             
%             p_no_drought = plot(period, data_GPP_no_drought, '-','Color', [9/255, 91/255, 208/255], 'LineWidth', 1.5);

            p_original = plot(period, data_GPP, '-k','LineWidth', 2);
            p_drought = plot(drought_period, drought_GPP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5); % 底部边界线
            
            p_no_drought = plot(period, data_GPP_no_drought, '-','Color',[9/255, 91/255, 208/255], 'LineWidth', 2.5);
            

            xlim_custom = [0 26];        
            xticks_custom = 0:1:26;      
            ylim_custom = [-2 14];    
            yticks_custom = -2:2:14;

            % 设置 x 和 y 轴范围
            xlim(xlim_custom);
            ylim(ylim_custom);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom);
            yticks(yticks_custom);

            set(gca, 'YTick', [0, 3, 6, 9, 12]);
            
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
            
            % 绘制新的时间点区间的蒸腾变化量曲线及其填充
            new_period_lush = [period(1:idx_start - 1), period_points(1), ceil(period_points(1)):floor(period_points(2)), period_points(2), period(idx_end:end)];
            new_GPP_rate_lush = [data_GPP(1:idx_start - 1), Interpolation_values_drought(1), interp1(period_points, Interpolation_values_drought, ceil(period_points(1)):floor(period_points(2))), Interpolation_values_drought(2), data_GPP(idx_end:end)];
            
            later_period = [period_points(2), period(idx_end:end)];
            later_GPP_rate = [Interpolation_values_drought(2), data_GPP(idx_end:end)];
            
            xverts_drought = [drought_period(1:end-1); drought_period(1:end-1); drought_period(2:end); drought_period(2:end)];
            yverts_drought = [zeros(1,length(drought_period)-1); drought_GPP_rate(1:end-1); drought_GPP_rate(2:end); zeros(1,length(drought_period)-1)];
            p_later = patch(xverts_drought, yverts_drought, [9/255, 91/255, 208/255], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
            
            % 构造填充区域的顶点坐标
            xverts_fill = [drought_period, fliplr(lush_period)];
            yverts_fill = [drought_GPP_rate, fliplr(lush_GPP_rate)];
            
            % 绘制顶部边界线
            hold on;

            % 计算光带的上下边界
            offset = (0.08*18)/16; % 定义偏移量，可以调整这个值来改变光带的宽度
            upper_bound = data_GPP_no_drought + offset;
            lower_bound = data_GPP_no_drought - offset;

            % 第一层光带，透明度较高
%             fill([period, fliplr(period)], [upper_bound+(0.08*18)/16, fliplr(lower_bound-(0.08*18)/16)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
%             % 第二层光带，透明度较低
%             fill([period, fliplr(period)], [upper_bound, fliplr(lower_bound)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

            % 绘制不干旱的情况
%             plot(period, data_GPP_no_drought, '-','Color',[12/255, 114/255, 31/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');
            plot(period, data_GPP_no_drought, '-','Color', [9/255, 91/255, 208/255], 'LineWidth', 2.5 ,'HandleVisibility', 'off');

            % 使用亮眼的颜色进行填充，并设置透明度
            % highlight_color = [184/255, 160/255, 52/255];%枯萎！！
            highlight_color = [66/255, 66/255, 148/255];
            p_fill = patch(xverts_fill, yverts_fill, highlight_color, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 无边界颜色

            plot(drought_period, drought_GPP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off'); % 底部边界线

            plot(drought_period([1, end]), drought_GPP_rate([1, end]), 'p', ...
                'MarkerSize', 3, 'LineWidth', 2, 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255], 'HandleVisibility', 'off');
            
            % 设置标题和 y 轴标签，根据 third_word 和 drought_speed
            if isempty(tree_kind)
                if strcmp(third_word, 'FD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'GPP')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = 'g C m⁻² d⁻¹';
                    end
                elseif strcmp(third_word, 'SD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Slow Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'GPP')
                        title_str = sprintf('%s Response to %s Slow Drought', data_kind, season_file);
                        ylabel_str = 'g C m⁻² d⁻¹';
                    end
                end
            else
                if strcmp(third_word, 'FD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'GPP')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'g C m⁻² d⁻¹';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'GPP')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'g C m⁻² d⁻¹';
                        end
                    end
                elseif strcmp(third_word, 'SD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'GPP')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'g C m⁻² d⁻¹';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'GPP')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'g C m⁻² d⁻¹';
                        end
                    end
                end
            end

            % 获取当前坐标轴对象
            ax = gca;
            set(gca,  'XTickLabel', []);
            
            % 获取坐标轴位置信息
            ax_position = ax.Position;
            
            ax.YAxis.Label.Position(1) = -1.5;
            
            % 设置 x 轴和 y 轴的线条粗细
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            ax.YAxis.FontWeight = 'bold';

                                                % 获取当前子图的当前位置
            pos = get(gca, 'Position');
            
            % 增加高度并保持上端位置不变
            new_height = pos(4) * 1.3;  % 增加高度
            pos(2) = pos(2) - (new_height - pos(4));  % 调整下端位置
            pos(4) = new_height;  % 设置新的高度
            
            % 应用新的位置
            set(gca, 'Position', pos);

            title_handle_1 = title(title_str, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            ylabel(ylabel_str, 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');


            if strcmp(season_file, 'Spring')
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Productivity', 'Productivity Anomaly'}, 'Location', 'northwest','FontSize', 14, 'FontWeight', 'bold');
            else
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Productivity', 'Productivity Anomaly'}, 'Location', 'northeast','FontSize', 18, 'FontWeight', 'bold');
            end

            legend_handle.Box = 'off';

            y_norm2 = 0.06/3;
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
            y_norm3 = 0.15/3;
            y_relative3 = y_norm3 * ax_position(4);

            box off     % 取消边框
            ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax1.XAxis.LineWidth = 2;
            ax1.YAxis.LineWidth = 2;
            ax1.YAxis.FontWeight = 'bold';

            

            % 第二个子图：绘制新的时间点区间的蒸腾变化量曲线和填充色
            h2 = subplot(2, 1, 2);
            plot(periods_before_drought, cumulative_transpiration_before_drought, '-k','LineWidth', 2.5);
            hold on
            
            offset = cumulative_transpiration_before_drought(end); % 你想加上的值
            
            % 修改后的数据
            cumulative_transpiration_drought_offset = cumulative_transpiration_droughting + offset;
            cumulative_transpiration_lush_offset = cumulative_transpiration_lush + offset;
            
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5);
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [9/255, 91/255, 208/255, 0.8],'LineWidth',2.5);

            % 根据 fourth_word 设置不同的 x 和 y 轴范围和刻度间隔
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 120];    
                            yticks_custom2 = 0:30:120;
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 120];    
                            yticks_custom2 = 0:30:120;
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        xlim_custom2 = [0 26];        
                        xticks_custom2 = 0:1:26;      
                        ylim_custom2 = [0 120];    
                        yticks_custom2 = 0:30:120;
                end
            end

            % 设置 x 和 y 轴范围
            xlim(xlim_custom2);
            ylim(ylim_custom2);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom2);
            yticks(yticks_custom2);
            
            
            % 创建阴影区域
            x_fill_shadow_1 = [period_points(1), serious_period, serious_period, period_points(1)];
            y_fill_shadow_1 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_1, y_fill_shadow_1, [128/255, 128/255, 128/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_2 = [serious_period, period_points(2), period_points(2), serious_period];
            y_fill_shadow_2 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_2, y_fill_shadow_2, [128/255, 128/255, 128/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_3 = [period_points(2), 26, 26, period_points(2)];
            y_fill_shadow_3 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_3, y_fill_shadow_3, [128/255, 128/255, 128/255], 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            hold on
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [9/255, 91/255, 208/255],'LineWidth',2.5, 'HandleVisibility', 'off');
            
            % 添加箭头、短线
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [9/255, 91/255, 208/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [9/255, 91/255, 208/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-1.5, '^r', ...
                            'MarkerFaceColor', [9/255, 91/255, 208/255], 'MarkerEdgeColor', [9/255, 91/255, 208/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+3, 'vr', ...
                            'MarkerFaceColor', [9/255, 91/255, 208/255], 'MarkerEdgeColor', [9/255, 91/255, 208/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [9/255, 91/255, 208/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-3, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+1.5, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);

                    end

                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [9/255, 91/255, 208/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [9/255, 91/255, 208/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-1.5, '^r', ...
                            'MarkerFaceColor', [9/255, 91/255, 208/255], 'MarkerEdgeColor', [9/255, 91/255, 208/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+3, 'vr', ...
                            'MarkerFaceColor', [9/255, 91/255, 208/255], 'MarkerEdgeColor', [9/255, 91/255, 208/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [9/255, 91/255, 208/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-3, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+1.5, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);

                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    % 对于不干旱的情况
                    % 绘制横线
                    line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [9/255, 91/255, 208/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [9/255, 91/255, 208/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(lush_period(1), cumulative_transpiration_lush_offset(end)-1.5, '^r', ...
                        'MarkerFaceColor', [9/255, 91/255, 208/255], 'MarkerEdgeColor', [9/255, 91/255, 208/255], 'MarkerSize', 5);  % 上箭头
                    plot(lush_period(1), cumulative_transpiration_lush_offset(1)+3, 'vr', ...
                        'MarkerFaceColor', [9/255, 91/255, 208/255], 'MarkerEdgeColor', [9/255, 91/255, 208/255],'MarkerSize', 5);  % 下箭头
                    
                    % 计算线段长度
                    length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                        ['{\bf ', num2str(round(length_of_line, 2)), '} \bf{}'], ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [9/255, 91/255, 208/255]);
                    
                    % 对于干旱情况
                    % 绘制横线
                    line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(drought_period(end), cumulative_transpiration_drought_offset(end)-3, '^r', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                    plot(drought_period(end), cumulative_transpiration_drought_offset(1)+2, 'vr', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                    
                    
                    % 计算线段长度
                    length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)])-0.05, ...
                        ['{\bf}', num2str(round(length_of_line2, 2))], ...
                        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);

                    
                end
            end
            
            scatter(drought_period(1), cumulative_transpiration_drought_offset(1), ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(lush_period(end), cumulative_transpiration_lush_offset(end), ...
                'MarkerFaceColor', [9/255, 91/255, 208/255], 'MarkerEdgeColor', [9/255, 91/255, 208/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            % grid on;
            % 黑色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+1, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');

                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+1, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        text(drought_period(1), cumulative_transpiration_drought_offset(1)+1, ...
                            ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                            'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                end
            end

            % 红色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.03, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.02, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.02, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                end
            end

            % 绿色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(lush_period(end), cumulative_transpiration_lush_offset(end)+5, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [9/255, 91/255, 208/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        text(lush_period(end), cumulative_transpiration_lush_offset(end)+5, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [9/255, 91/255, 208/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                     text(lush_period(end), cumulative_transpiration_lush_offset(end)+6, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [9/255, 91/255, 208/255]);
                end
            end

            
            % 获取当前坐标轴对象
            ax2 = gca;
            set(gca,  'XTickLabel', []);

            % 设置 x 轴和 y 轴的线条粗细
            ax2.XAxis.LineWidth = 2;
            ax2.YAxis.LineWidth = 2;
            ax2.XAxis.FontSize = 20;
            ax2.YAxis.FontSize = 20;
            ax2.YAxis.FontWeight = 'bold';

            title_handle_2 = title('Cumulative Gross Primary Productivity', 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
            ylabel('g C m⁻²', 'FontName', 'Times New Roman','FontSize', 24, 'FontWeight', 'bold');
            
            if strcmp(season_file, 'Spring')
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 14, 'FontWeight', 'bold');
            else
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northeast', 'FontSize', 18, 'FontWeight', 'bold');
            end
            h_legend.Box = 'off';
            
            % 获取坐标轴位置信息
            ax2_position2 = ax2.Position;

            % 调整 x 轴标签位置
            y_norm2 = 0.075/3 ;
            y_relative = y_norm2 * ax_position(4);

            abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            % 在原有刻度上添加额外的标记
            for i = 1:length(positions)
                x_position = positions(i);
                x_norm2 = (x_position - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));

                x_relative = x_norm2 * ax2_position2(3);

                % 绘制垂直于 x 轴的线条
                line([x_relative*(26/ax2_position2(3)), x_relative*(26/ax2_position2(3))], [abs_length + ylim_custom2(1),ylim_custom2(1)], ...
                    'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            
            % 调整 x 轴标签位置
            abs_length3 = y_relative3*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            ax2.XAxis.Label.Position(2) = ax2.XAxis.Label.Position(2) - abs_length3;

            % 循环添加标签
            for i = 1:length(positions)
                x_label = positions(i);
                x_norm_label = (x_label - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));
        
                x_relative_label = x_norm_label * ax2_position2(3);
        
                abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));
                
                % 添加标签
                text(x_relative_label*(26/ax2_position2(3)), ylim_custom2(1)-abs_length, labels{i}, 'Color', 'k', 'VerticalAlignment', 'top', ...
                    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
            end
            
            ax2.YAxis.Label.Position(1) = -1.5;
            
            % 调整垂直间距
            vertical_gap = 0.06;
            % 获取第一个子图的位置
            pos1 = get(h1, 'Position');
            
            % 获取第二个子图的位置，并调整其位置
            pos2 = get(h2, 'Position');
            pos2(2) = pos1(2) - pos1(4) + vertical_gap; % 调整第二个子图的 y 坐标位置
            
            % 更新第二个子图的位置
            set(h2, 'Position', pos2);

            box off     % 取消边框
            ax2_2 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax2_2,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax2_2.XAxis.LineWidth = 2;
            ax2_2.YAxis.LineWidth = 2;
            ax2_2.YAxis.FontWeight = 'bold';
            
            % 设置保存路径和文件名
            save_path = ['E:\Data\Drought_quantification\Quantify(pictures)\Quantify_0403\' segmented_name];
            
            % 完整的保存路径
            full_file_path = fullfile(save_path, output_file);
            
            % 设置输出分辨率为 300 DPI
            resolution = 500;
            
            % 保存图形为 JPEG 格式
            saveas(gcf, full_file_path, 'svg');

        elseif strcmp(data_kind, 'NEP')
            data_NEP = readmatrix(selected_file);
            data_NEP = data_NEP(:)';

            % 获取 period_points
            for n = 1:numel(selected_files_periods)
                selected_file_period = selected_files_periods{n};
                [~, period_kinds, ~] = fileparts(selected_file_period);
                period_words = split(period_kinds, '_');
                period_kind = period_words{4};
    
                if strcmp(period_kind, 'NEP') % 假设period_kind也应该是NEP
                    period_points = readmatrix(selected_file_period);
                    period_points = period_points(:)';
                    break; % 找到匹配的文件就退出循环
                end
            end

            Interpolation_values_drought = interp1(period, data_NEP, period_points);

            selected_file_path_parts = split(selected_file, filesep);
            fourth_file_name = selected_file_path_parts(5);

            if strcmp(fourth_file_name, 'Season_ALL_1_month')
                % 构造对应文件的路径
                corresponding_file = strrep(selected_file, 'Season_ALL_1_month', 'Season_ALL_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(10) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s.svg', data_kind, third_word, season_file);
            else
                corresponding_file = strrep(selected_file, 'Season_N&P_1_month', 'Season_N&P_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(11) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s(%s).svg', data_kind, third_word, season_file, tree_kind);
            end
            
            % 重新组合路径
            corresponding_file = fullfile(path_parts{:});

            data_NEP_no_drought = readmatrix(corresponding_file);
            data_NEP_no_drought = data_NEP_no_drought(:)';

            Interpolation_values_no_drought = interp1(period, data_NEP_no_drought, period_points);

            % 提取从 start_period 开始的时间和蒸腾变化量数据
            idx_start = find(period >= ceil(period_points(1)), 1);
            idx_end = find(period >= ceil(period_points(2)), 1);
            drought_period = [period_points(1), period(idx_start:idx_end-1),period_points(2)];
            drought_NEP_rate = [Interpolation_values_drought(1), data_NEP(idx_start:idx_end-1),Interpolation_values_drought(2)];
            
            % 干旱期间的累计蒸腾量
            cumulative_transpiration_droughting = cumtrapz(drought_period, drought_NEP_rate);
            
            new_drought_period = [drought_period, period(idx_end:end)];
            new_drought_NEP_rate = [drought_NEP_rate, data_NEP(idx_end:end)];
            
            cumulative_transpiration_drought = cumtrapz(new_drought_period, new_drought_NEP_rate);
            
            periods_before_drought = [period(1:idx_start - 1), period_points(1)];
            NEP_rate_before_drought = [data_NEP(1:idx_start - 1), Interpolation_values_drought(1)];
            cumulative_transpiration_before_drought = cumtrapz(periods_before_drought, NEP_rate_before_drought);
            
            new_period = [period(1:idx_start - 1), drought_period, period(idx_end:end)];
            new_NEP_rate = [data_NEP(1:idx_start - 1), drought_NEP_rate, data_NEP(idx_end:end)];
            
            % 假设没有发生干旱 假设没有发生干旱 假设没有发生干旱
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况

            integer_periods = ceil(period_points(1)):floor(period_points(2));
            
            % 计算整个时间范围内的累积蒸腾变化量
            cumulative_transpiration_total = cumtrapz(new_period, new_NEP_rate);
            
            lush_period = [period_points(1), integer_periods, period_points(2)];
            lush_NEP_rate = [Interpolation_values_no_drought(1), data_NEP_no_drought(integer_periods), Interpolation_values_no_drought(2)]; 
            
            % 计算没有干旱发生 的累计蒸腾量 针对干旱时期对应的时间
            cumulative_transpiration_lush = cumtrapz(lush_period, lush_NEP_rate);
            
            new_lush_period = [lush_period, period(idx_end:end)];
            new_lush_NEP_rate = [lush_NEP_rate, data_NEP_no_drought(idx_end:end)];
           
            % 计算没有干旱发生 的累计蒸腾量 针对干旱开始到25！！
            cumulative_transpiration_lush_after_drought = cumtrapz(new_lush_period, new_lush_NEP_rate);
            
            new_lush_period_total = [period(1:idx_start - 1), new_lush_period];
            new_lush_NEP_rate_total = [data_NEP_no_drought(1:idx_start - 1), new_lush_NEP_rate];
            
            cumulative_transpiration_lush_total = cumtrapz(new_lush_period_total, new_lush_NEP_rate_total);
            
            % 创建图形窗口并调整尺寸
            % figure('Units', 'normalized', 'Position', [0.6, 0.1, 0.3, 0.5]);
            figure('Units', 'pixels', 'Position', [1750,100,800,1000]);
            % 绘制蒸腾变化量随时间变化的曲线及其填充
            h1 = subplot(2, 1, 1);
            % 绘制干旱的情况
            p_star = plot(drought_period([1, end]), drought_NEP_rate([1, end]), 'p', 'MarkerSize', 3, 'LineWidth', 2 ,...
                'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255]);
            hold on;
%             p_original = plot(period, data_NEP, '-','Color', [23/255, 216/255, 15/255, 0.5],'LineWidth', 1.5);
%             
%             
%             p_drought = plot(drought_period, drought_NEP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.5); % 底部边界线
%             
% 
%             p_no_drought = plot(period, data_NEP_no_drought, '-','Color',[23/255, 216/255, 15/255], 'LineWidth', 1.5);

            p_original = plot(period, data_NEP, '-k','LineWidth', 2);
            p_drought = plot(drought_period, drought_NEP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5); % 底部边界线
            
            p_no_drought = plot(period, data_NEP_no_drought, '-','Color',[23/255, 216/255, 15/255], 'LineWidth', 2.5);
            

            xlim_custom = [0 26];        
            xticks_custom = 0:1:26;      
            ylim_custom = [-1.5 2];
            yticks_custom = -1.5:0.5:2;

            % 设置 x 和 y 轴范围
            xlim(xlim_custom);
            ylim(ylim_custom);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom);
            yticks(yticks_custom);

            set(gca, 'YTick', [-1.0, -0.5, 0, 0.5, 1.0, 1.5]);

            
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
            
            % 绘制新的时间点区间的蒸腾变化量曲线及其填充
            new_period_lush = [period(1:idx_start - 1), period_points(1), ceil(period_points(1)):floor(period_points(2)), period_points(2), period(idx_end:end)];
            new_NEP_rate_lush = [data_NEP(1:idx_start - 1), Interpolation_values_drought(1), interp1(period_points, Interpolation_values_drought, ceil(period_points(1)):floor(period_points(2))), Interpolation_values_drought(2), data_NEP(idx_end:end)];
            
            later_period = [period_points(2), period(idx_end:end)];
            later_NEP_rate = [Interpolation_values_drought(2), data_NEP(idx_end:end)];
            
            xverts_drought = [drought_period(1:end-1); drought_period(1:end-1); drought_period(2:end); drought_period(2:end)];
            yverts_drought = [zeros(1,length(drought_period)-1); drought_NEP_rate(1:end-1); drought_NEP_rate(2:end); zeros(1,length(drought_period)-1)];
%             p_later = patch(xverts_drought, yverts_drought, [23/255, 216/255, 15/255], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
            
            % 构造填充区域的顶点坐标
            xverts_fill = [drought_period, fliplr(lush_period)];
            yverts_fill = [drought_NEP_rate, fliplr(lush_NEP_rate)];
            
            % 绘制底部边界线
            hold on;

            % 计算光带的上下边界
            offset = (0.08*3)/16; % 定义偏移量，可以调整这个值来改变光带的宽度
            upper_bound = data_NEP_no_drought + offset;
            lower_bound = data_NEP_no_drought - offset;

            % 第一层光带，透明度较高
%             fill([period, fliplr(period)], [upper_bound+(0.08*3)/16, fliplr(lower_bound-(0.08*3)/16)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
%             % 第二层光带，透明度较低
%             fill([period, fliplr(period)], [upper_bound, fliplr(lower_bound)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

            % 绘制不干旱的情况
            plot(period, data_NEP_no_drought, '-','Color',[23/255, 216/255, 15/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');

            % 使用亮眼的颜色进行填充，并设置透明度
            % highlight_color = [184/255, 160/255, 52/255];%枯萎！！
            highlight_color = [66/255, 66/255, 148/255];
            p_fill = patch(xverts_fill, yverts_fill, highlight_color, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 无边界颜色

            plot(drought_period, drought_NEP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off'); % 底部边界线

            plot(drought_period([1, end]), drought_NEP_rate([1, end]), 'p', ...
                'MarkerSize', 3, 'LineWidth', 2, 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255], 'HandleVisibility', 'off');
            
            % 设置标题和 y 轴标签，根据 third_word 和 drought_speed
            if isempty(tree_kind)
                if strcmp(third_word, 'FD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'NEP')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = '';
                    end
                elseif strcmp(third_word, 'SD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Slow Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'NEP')
                        title_str = sprintf('%s Response to %s Slow Drought', data_kind, season_file);
                        ylabel_str = '';
                    end
                end
            else
                if strcmp(third_word, 'FD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NEP')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NEP')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    end
                elseif strcmp(third_word, 'SD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NEP')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NEP')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    end
                end
            end

            % 获取当前坐标轴对象
            ax = gca;
            set(gca,  'XTickLabel', []);
            
            % 获取坐标轴位置信息
            ax_position = ax.Position;
            
            ax.YAxis.Label.Position(1) = -1.5;
            
            % 设置 x 轴和 y 轴的线条粗细
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            ax.YAxis.FontWeight = 'bold';

                                                % 获取当前子图的当前位置
            pos = get(gca, 'Position');
            
            % 增加高度并保持上端位置不变
            new_height = pos(4) * 1.3;  % 增加高度
            pos(2) = pos(2) - (new_height - pos(4));  % 调整下端位置
            pos(4) = new_height;  % 设置新的高度
            
            % 应用新的位置
            set(gca, 'Position', pos);

            title_handle_1 = title(title_str, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            ylabel(ylabel_str, 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
            if strcmp(season_file, 'Spring') || strcmp(season_file, 'Autumn')
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Productivity Anomaly'}, 'Location', 'northwest','FontSize', 18, 'FontWeight', 'bold');
            else
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Productivity Anomaly'}, 'Location', 'northeast','FontSize', 14, 'FontWeight', 'bold');
            end
            legend_handle.Box = 'off';

            y_norm2 = 0.06/3;
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
            y_norm3 = 0.15/3;
            y_relative3 = y_norm3 * ax_position(4);

            box off     % 取消边框
            ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax1.XAxis.LineWidth = 2;
            ax1.YAxis.LineWidth = 2;
            ax1.YAxis.FontWeight = 'bold';

            

            % 第二个子图：绘制新的时间点区间的蒸腾变化量曲线和填充色
            h2 = subplot(2, 1, 2);
            plot(periods_before_drought, cumulative_transpiration_before_drought, '-k','LineWidth', 2.5);
            hold on
            
            offset = cumulative_transpiration_before_drought(end); % 你想加上的值
            
            % 修改后的数据
            cumulative_transpiration_drought_offset = cumulative_transpiration_droughting + offset;
            cumulative_transpiration_lush_offset = cumulative_transpiration_lush + offset;
            
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5);
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [23/255, 216/255, 15/255],'LineWidth',2.5);

            % 根据 fourth_word 设置不同的 x 和 y 轴范围和刻度间隔
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [-3 1];    
                            yticks_custom2 = -3:1:1;
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [-5 1];    
                            yticks_custom2 = -5:2:1;
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        xlim_custom2 = [0 26];        
                        xticks_custom2 = 0:1:26;      
                        ylim_custom2 = [-3 1];    
                        yticks_custom2 = -3:1:1;
                end
            end

            % 设置 x 和 y 轴范围
            xlim(xlim_custom2);
            ylim(ylim_custom2);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom2);
            yticks(yticks_custom2);
            
            % 创建阴影区域
            x_fill_shadow_1 = [period_points(1), serious_period, serious_period, period_points(1)];
            y_fill_shadow_1 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_1, y_fill_shadow_1, [128/255, 128/255, 128/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_2 = [serious_period, period_points(2), period_points(2), serious_period];
            y_fill_shadow_2 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_2, y_fill_shadow_2, [128/255, 128/255, 128/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_3 = [period_points(2), 26, 26, period_points(2)];
            y_fill_shadow_3 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_3, y_fill_shadow_3, [128/255, 128/255, 128/255], 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            hold on
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [23/255, 216/255, 15/255],'LineWidth',2.5, 'HandleVisibility', 'off');
            
            % 添加箭头、短线
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        
                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [23/255, 216/255, 15/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [23/255, 216/255, 15/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.05, '^r', ...
                            'MarkerFaceColor', [23/255, 216/255, 15/255], 'MarkerEdgeColor', [23/255, 216/255, 15/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+0.1, 'vr', ...
                            'MarkerFaceColor', [23/255, 216/255, 15/255], 'MarkerEdgeColor', [23/255, 216/255, 15/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [23/255, 216/255, 15/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-0.35, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.4, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);

                       
                    end

                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        

                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [23/255, 216/255, 15/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [23/255, 216/255, 15/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.07, '^r', ...
                            'MarkerFaceColor', [23/255, 216/255, 15/255], 'MarkerEdgeColor', [23/255, 216/255, 15/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+0.1, 'vr', ...
                            'MarkerFaceColor', [23/255, 216/255, 15/255], 'MarkerEdgeColor', [23/255, 216/255, 15/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)])+0.05, ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [23/255, 216/255, 15/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-0.1, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.08, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)])+0.05, ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    

                    case 'Autumn'
                    % 对于不干旱的情况
                    % 绘制横线
                    line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [23/255, 216/255, 15/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [23/255, 216/255, 15/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.05, '^r', ...
                        'MarkerFaceColor', [23/255, 216/255, 15/255], 'MarkerEdgeColor', [23/255, 216/255, 15/255], 'MarkerSize', 5);  % 上箭头
                    plot(lush_period(1), cumulative_transpiration_lush_offset(1)+0.1, 'vr', ...
                        'MarkerFaceColor', [23/255, 216/255, 15/255], 'MarkerEdgeColor', [23/255, 216/255, 15/255],'MarkerSize', 5);  % 下箭头
                    
                    % 计算线段长度
                    length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                        ['{\bf ', num2str(round(length_of_line, 2)), '} \bf{}'], ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [23/255, 216/255, 15/255]);
                    
                    % 对于干旱情况
                    % 绘制横线
                    line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(drought_period(end), cumulative_transpiration_drought_offset(end)-0.1, '^r', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                    plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.1, 'vr', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                    
                    
                    % 计算线段长度
                    length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(drought_period(end)+0.15, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)])-0.05, ...
                        ['{\bf—————}', num2str(round(length_of_line2, 2))], ...
                        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);

                    
                end
            end
            
            scatter(drought_period(1), cumulative_transpiration_drought_offset(1), ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(lush_period(end), cumulative_transpiration_lush_offset(end), ...
                'MarkerFaceColor', [23/255, 216/255, 15/255], 'MarkerEdgeColor', [23/255, 216/255, 15/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            % grid on;
            % 黑色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                       
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.035, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.1, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.04, ...
                            ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                            'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                end
            end

            % 红色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.035, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.02, ...
                            ['\bf{——} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.02, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                end
            end

            % 绿色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(lush_period(end), cumulative_transpiration_lush_offset(end)+0.035, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [23/255, 216/255, 15/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        text(lush_period(end), cumulative_transpiration_lush_offset(end)+0.03, ...
                            ['\bf{——} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [23/255, 216/255, 15/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                     text(lush_period(end), cumulative_transpiration_lush_offset(end)+0.03, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [23/255, 216/255, 15/255]);
                end
            end

            
            % 获取当前坐标轴对象
            ax2 = gca;
            set(gca,  'XTickLabel', []);

            % 设置 x 轴和 y 轴的线条粗细
            ax2.XAxis.LineWidth = 2;
            ax2.YAxis.LineWidth = 2;
            ax2.XAxis.FontSize = 20;
            ax2.YAxis.FontSize = 20;
            ax2.YAxis.FontWeight = 'bold';

            title_handle_2 = title('Cumulative Net Ecosystem Productivity', 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
            ylabel('', 'FontName', 'Times New Roman','FontSize', 24, 'FontWeight', 'bold');
            
            if strcmp(season_file, 'Spring')
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 14, 'FontWeight', 'bold');
            else
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northeast', 'FontSize', 18, 'FontWeight', 'bold');
            end
            h_legend.Box = 'off';
            
            % 获取坐标轴位置信息
            ax2_position2 = ax2.Position;

                                    % 调整 x 轴标签位置
            y_norm2 = 0.075/3 ;
            y_relative = y_norm2 * ax_position(4);

            abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            % 在原有刻度上添加额外的标记
            for i = 1:length(positions)
                x_position = positions(i);
                x_norm2 = (x_position - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));

                x_relative = x_norm2 * ax2_position2(3);

                % 绘制垂直于 x 轴的线条
                line([x_relative*(26/ax2_position2(3)), x_relative*(26/ax2_position2(3))], [abs_length + ylim_custom2(1),ylim_custom2(1)], ...
                    'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            
            % 调整 x 轴标签位置
            abs_length3 = y_relative3*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            ax2.XAxis.Label.Position(2) = ax2.XAxis.Label.Position(2) - abs_length3;

            % 循环添加标签
            for i = 1:length(positions)
                x_label = positions(i);
                x_norm_label = (x_label - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));
        
                x_relative_label = x_norm_label * ax2_position2(3);
        
                abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));
                
                % 添加标签
                text(x_relative_label*(26/ax2_position2(3)), ylim_custom2(1)-abs_length, labels{i}, 'Color', 'k', 'VerticalAlignment', 'top', ...
                    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
            end
            
            ax2.YAxis.Label.Position(1) = -1.5;
            
            % 调整垂直间距
            vertical_gap = 0.06;
            % 获取第一个子图的位置
            pos1 = get(h1, 'Position');
            
            % 获取第二个子图的位置，并调整其位置
            pos2 = get(h2, 'Position');
            pos2(2) = pos1(2) - pos1(4) + vertical_gap; % 调整第二个子图的 y 坐标位置
            
            % 更新第二个子图的位置
            set(h2, 'Position', pos2);

            box off     % 取消边框
            ax2_2 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax2_2,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax2_2.XAxis.LineWidth = 2;
            ax2_2.YAxis.LineWidth = 2;
            ax2_2.YAxis.FontWeight = 'bold';
            
            % 设置保存路径和文件名
            save_path = ['E:\Data\Drought_quantification\Quantify(pictures)\Quantify_0403\' segmented_name];
            
            % 完整的保存路径
            full_file_path = fullfile(save_path, output_file);
            
            % 设置输出分辨率为 300 DPI
            resolution = 500;
            
            % 保存图形为 JPEG 格式
            saveas(gcf, full_file_path, 'svg');

        elseif strcmp(data_kind, 'NPP')
            data_NPP = readmatrix(selected_file);
            data_NPP = data_NPP(:)';

            % 获取 period_points
            for n = 1:numel(selected_files_periods)
                selected_file_period = selected_files_periods{n};
                [~, period_kinds, ~] = fileparts(selected_file_period);
                period_words = split(period_kinds, '_');
                period_kind = period_words{4};
    
                if strcmp(period_kind, 'NPP') % 假设period_kind也应该是NPP
                    period_points = readmatrix(selected_file_period);
                    period_points = period_points(:)';
                    break; % 找到匹配的文件就退出循环
                end
            end

            Interpolation_values_drought = interp1(period, data_NPP, period_points);

            selected_file_path_parts = split(selected_file, filesep);
            fourth_file_name = selected_file_path_parts(5);

            if strcmp(fourth_file_name, 'Season_ALL_1_month')
                % 构造对应文件的路径
                corresponding_file = strrep(selected_file, 'Season_ALL_1_month', 'Season_ALL_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(10) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s.svg', data_kind, third_word, season_file);
            else
                corresponding_file = strrep(selected_file, 'Season_N&P_1_month', 'Season_N&P_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(11) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s(%s).svg', data_kind, third_word, season_file, tree_kind);
            end
            
            % 重新组合路径
            corresponding_file = fullfile(path_parts{:});

            data_NPP_no_drought = readmatrix(corresponding_file);
            data_NPP_no_drought = data_NPP_no_drought(:)';

            Interpolation_values_no_drought = interp1(period, data_NPP_no_drought, period_points);

            % 提取从 start_period 开始的时间和蒸腾变化量数据
            idx_start = find(period >= ceil(period_points(1)), 1);
            idx_end = find(period >= ceil(period_points(2)), 1);
            drought_period = [period_points(1), period(idx_start:idx_end-1),period_points(2)];
            drought_NPP_rate = [Interpolation_values_drought(1), data_NPP(idx_start:idx_end-1),Interpolation_values_drought(2)];
            
            % 干旱期间的累计蒸腾量
            cumulative_transpiration_droughting = cumtrapz(drought_period, drought_NPP_rate);
            
            new_drought_period = [drought_period, period(idx_end:end)];
            new_drought_NPP_rate = [drought_NPP_rate, data_NPP(idx_end:end)];
            
            cumulative_transpiration_drought = cumtrapz(new_drought_period, new_drought_NPP_rate);
            
            periods_before_drought = [period(1:idx_start - 1), period_points(1)];
            NPP_rate_before_drought = [data_NPP(1:idx_start - 1), Interpolation_values_drought(1)];
            cumulative_transpiration_before_drought = cumtrapz(periods_before_drought, NPP_rate_before_drought);
            
            new_period = [period(1:idx_start - 1), drought_period, period(idx_end:end)];
            new_NPP_rate = [data_NPP(1:idx_start - 1), drought_NPP_rate, data_NPP(idx_end:end)];
            
            % 假设没有发生干旱 假设没有发生干旱 假设没有发生干旱
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况

            integer_periods = ceil(period_points(1)):floor(period_points(2));
            
            % 计算整个时间范围内的累积蒸腾变化量
            cumulative_transpiration_total = cumtrapz(new_period, new_NPP_rate);
            
            lush_period = [period_points(1), integer_periods, period_points(2)];
            lush_NPP_rate = [Interpolation_values_no_drought(1), data_NPP_no_drought(integer_periods), Interpolation_values_no_drought(2)]; 
            
            % 计算没有干旱发生 的累计蒸腾量 针对干旱时期对应的时间
            cumulative_transpiration_lush = cumtrapz(lush_period, lush_NPP_rate);
            
            new_lush_period = [lush_period, period(idx_end:end)];
            new_lush_NPP_rate = [lush_NPP_rate, data_NPP_no_drought(idx_end:end)];
           
            % 计算没有干旱发生 的累计蒸腾量 针对干旱开始到25！！
            cumulative_transpiration_lush_after_drought = cumtrapz(new_lush_period, new_lush_NPP_rate);
            
            new_lush_period_total = [period(1:idx_start - 1), new_lush_period];
            new_lush_NPP_rate_total = [data_NPP_no_drought(1:idx_start - 1), new_lush_NPP_rate];
            
            cumulative_transpiration_lush_total = cumtrapz(new_lush_period_total, new_lush_NPP_rate_total);
            
            % 创建图形窗口并调整尺寸
            % figure('Units', 'normalized', 'Position', [0.6, 0.1, 0.3, 0.5]);
            figure('Units', 'pixels', 'Position', [1750,100,800,1000]);
            % 绘制蒸腾变化量随时间变化的曲线及其填充
            h1 = subplot(2, 1, 1);
            % 绘制干旱的情况
            p_star = plot(drought_period([1, end]), drought_NPP_rate([1, end]), 'p', 'MarkerSize', 3, 'LineWidth', 2 ,...
                'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255]);
            hold on;

%             p_original = plot(period, data_NPP, '-','Color', [154/255, 104/255, 196/255, 0.8],'LineWidth', 1.5);
%   
%             p_drought = plot(drought_period, drought_NPP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.5); % 底部边界线
% 
%             
%             p_no_drought = plot(period, data_NPP_no_drought, '-','Color', [154/255, 104/255, 196/255], 'LineWidth', 1.5);

            p_original = plot(period, data_NPP, '-k','LineWidth', 2);
            p_drought = plot(drought_period, drought_NPP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5); % 底部边界线
            
            p_no_drought = plot(period, data_NPP_no_drought, '-','Color',[154/255, 104/255, 196/255], 'LineWidth', 2.5);

            xlim_custom = [0 26];        
            xticks_custom = 0:1:26;      
            ylim_custom = [-1 6];    
            yticks_custom = -1:2:6;

            % 设置 x 和 y 轴范围
            xlim(xlim_custom);
            ylim(ylim_custom);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom);
            yticks(yticks_custom);

            set(gca, 'YTick', [0, 1, 2, 3, 4, 5]);
            
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
            
            % 绘制新的时间点区间的蒸腾变化量曲线及其填充
            new_period_lush = [period(1:idx_start - 1), period_points(1), ceil(period_points(1)):floor(period_points(2)), period_points(2), period(idx_end:end)];
            new_NPP_rate_lush = [data_NPP(1:idx_start - 1), Interpolation_values_drought(1), interp1(period_points, Interpolation_values_drought, ceil(period_points(1)):floor(period_points(2))), Interpolation_values_drought(2), data_NPP(idx_end:end)];
            
            later_period = [period_points(2), period(idx_end:end)];
            later_NPP_rate = [Interpolation_values_drought(2), data_NPP(idx_end:end)];
            
            xverts_drought = [drought_period(1:end-1); drought_period(1:end-1); drought_period(2:end); drought_period(2:end)];
            yverts_drought = [zeros(1,length(drought_period)-1); drought_NPP_rate(1:end-1); drought_NPP_rate(2:end); zeros(1,length(drought_period)-1)];
            p_later = patch(xverts_drought, yverts_drought, [154/255, 104/255, 196/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            
            % 构造填充区域的顶点坐标
            xverts_fill = [drought_period, fliplr(lush_period)];
            yverts_fill = [drought_NPP_rate, fliplr(lush_NPP_rate)];
            
            % 绘制顶部边界线
            hold on;

            % 计算光带的上下边界
            offset = (0.08*7)/16; % 定义偏移量，可以调整这个值来改变光带的宽度
            upper_bound = data_NPP_no_drought + offset;
            lower_bound = data_NPP_no_drought - offset;

%             % 第一层光带，透明度较高
%             fill([period, fliplr(period)], [upper_bound+(0.08*7)/16, fliplr(lower_bound-(0.08*7)/16)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
%             % 第二层光带，透明度较低
%             fill([period, fliplr(period)], [upper_bound, fliplr(lower_bound)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

            % 绘制不干旱的情况
            plot(period, data_NPP_no_drought, '-','Color',[154/255, 104/255, 196/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');

            % 使用亮眼的颜色进行填充，并设置透明度
            % highlight_color = [184/255, 160/255, 52/255];%枯萎！！
            highlight_color = [66/255, 66/255, 148/255];
            p_fill = patch(xverts_fill, yverts_fill, highlight_color, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 无边界颜色

            plot(drought_period, drought_NPP_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off'); % 底部边界线

            plot(drought_period([1, end]), drought_NPP_rate([1, end]), 'p', ...
                'MarkerSize', 3, 'LineWidth', 2, 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255], 'HandleVisibility', 'off');
            
            % 设置标题和 y 轴标签，根据 third_word 和 drought_speed
            if isempty(tree_kind)
                if strcmp(third_word, 'FD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'NPP')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = '';
                    end
                end
            else
                if strcmp(third_word, 'FD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NPP')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NPP')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    end
                end
            end

            % 获取当前坐标轴对象
            ax = gca;
            set(gca,  'XTickLabel', []);
            
            % 获取坐标轴位置信息
            ax_position = ax.Position;
            
            ax.YAxis.Label.Position(1) = -1.5;
            
            % 设置 x 轴和 y 轴的线条粗细
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            ax.YAxis.FontWeight = 'bold';

                                                % 获取当前子图的当前位置
            pos = get(gca, 'Position');
            
            % 增加高度并保持上端位置不变
            new_height = pos(4) * 1.3;  % 增加高度
            pos(2) = pos(2) - (new_height - pos(4));  % 调整下端位置
            pos(4) = new_height;  % 设置新的高度
            
            % 应用新的位置
            set(gca, 'Position', pos);

            title_handle_1 = title(title_str, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            ylabel(ylabel_str, 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');

            if strcmp(season_file, 'Spring')
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Productivity', 'Productivity Anomaly'}, 'Location', 'northwest','FontSize', 14, 'FontWeight', 'bold');
            else
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Productivity', 'Productivity Anomaly'}, 'Location', 'northeast','FontSize', 18, 'FontWeight', 'bold');
            end
            legend_handle.Box = 'off';

            y_norm2 = 0.06/3;
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
            y_norm3 = 0.15/3;
            y_relative3 = y_norm3 * ax_position(4);

            box off     % 取消边框
            ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax1.XAxis.LineWidth = 2;
            ax1.YAxis.LineWidth = 2;
            ax1.YAxis.FontWeight = 'bold';

            

            % 第二个子图：绘制新的时间点区间的蒸腾变化量曲线和填充色
            h2 = subplot(2, 1, 2);
            plot(periods_before_drought, cumulative_transpiration_before_drought, '-k','LineWidth', 2.5);
            hold on
            
            offset = cumulative_transpiration_before_drought(end); % 你想加上的值
            
            % 修改后的数据
            cumulative_transpiration_drought_offset = cumulative_transpiration_droughting + offset;
            cumulative_transpiration_lush_offset = cumulative_transpiration_lush + offset;
            
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5);
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [154/255, 104/255, 196/255],'LineWidth',2.5);

            % 根据 fourth_word 设置不同的 x 和 y 轴范围和刻度间隔
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 45];    
                            yticks_custom2 = 0:15:45;
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 45];    
                            yticks_custom2 = 0:15:45;
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        xlim_custom2 = [0 26];        
                        xticks_custom2 = 0:1:26;      
                        ylim_custom2 = [0 45];    
                        yticks_custom2 = 0:15:45;
                end
            end


            % 设置 x 和 y 轴范围
            xlim(xlim_custom2);
            ylim(ylim_custom2);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom2);
            yticks(yticks_custom2);
            
            % 创建阴影区域
            x_fill_shadow_1 = [period_points(1), serious_period, serious_period, period_points(1)];
            y_fill_shadow_1 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_1, y_fill_shadow_1, [128/255, 128/255, 128/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_2 = [serious_period, period_points(2), period_points(2), serious_period];
            y_fill_shadow_2 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_2, y_fill_shadow_2, [128/255, 128/255, 128/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_3 = [period_points(2), 26, 26, period_points(2)];
            y_fill_shadow_3 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_3, y_fill_shadow_3, [128/255, 128/255, 128/255], 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            hold on
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [154/255, 104/255, 196/255],'LineWidth',2.5, 'HandleVisibility', 'off');
            
            % 添加箭头、短线
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file

                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [154/255, 104/255, 196/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [154/255, 104/255, 196/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.6, '^r', ...
                            'MarkerFaceColor', [154/255, 104/255, 196/255], 'MarkerEdgeColor', [154/255, 104/255, 196/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+1, 'vr', ...
                            'MarkerFaceColor', [154/255, 104/255, 196/255], 'MarkerEdgeColor', [154/255, 104/255, 196/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [154/255, 104/255, 196/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-1, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.6, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);

                    end

                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [154/255, 104/255, 196/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [154/255, 104/255, 196/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.7, '^r', ...
                            'MarkerFaceColor', [154/255, 104/255, 196/255], 'MarkerEdgeColor', [154/255, 104/255, 196/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+1, 'vr', ...
                            'MarkerFaceColor', [154/255, 104/255, 196/255], 'MarkerEdgeColor', [154/255, 104/255, 196/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [154/255, 104/255, 196/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-1.1, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.6, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    % 对于不干旱的情况
                    % 绘制横线
                    line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [154/255, 104/255, 196/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [154/255, 104/255, 196/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.7, '^r', ...
                        'MarkerFaceColor', [154/255, 104/255, 196/255], 'MarkerEdgeColor', [154/255, 104/255, 196/255], 'MarkerSize', 5);  % 上箭头
                    plot(lush_period(1), cumulative_transpiration_lush_offset(1)+1.2, 'vr', ...
                        'MarkerFaceColor', [154/255, 104/255, 196/255], 'MarkerEdgeColor', [154/255, 104/255, 196/255],'MarkerSize', 5);  % 下箭头
                    
                    % 计算线段长度
                    length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                        ['{\bf ', num2str(round(length_of_line, 2)), '} \bf{}'], ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [154/255, 104/255, 196/255]);
                    
                    % 对于干旱情况
                    % 绘制横线
                    line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(drought_period(end), cumulative_transpiration_drought_offset(end)-1.2, '^r', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                    plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.7, 'vr', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                    
                    
                    % 计算线段长度
                    length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                        ['{\bf}', num2str(round(length_of_line2, 2))], ...
                        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                end
            end
            
            scatter(drought_period(1), cumulative_transpiration_drought_offset(1), ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(lush_period(end), cumulative_transpiration_lush_offset(end), ...
                'MarkerFaceColor', [154/255, 104/255, 196/255], 'MarkerEdgeColor', [154/255, 104/255, 196/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            % grid on;
            % 黑色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.5, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.3, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.25, ...
                            ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                            'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                end
            end

            % 红色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.04, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.02, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                end
            end

            % 绿色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(lush_period(end), cumulative_transpiration_lush_offset(end)+1.5, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [154/255, 104/255, 196/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                         text(lush_period(end), cumulative_transpiration_lush_offset(end), ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [154/255, 104/255, 196/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                     text(lush_period(end), cumulative_transpiration_lush_offset(end)+1.5, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [154/255, 104/255, 196/255]);
                end
            end

            
            % 获取当前坐标轴对象
            ax2 = gca;
            set(gca,  'XTickLabel', []);

            % 设置 x 轴和 y 轴的线条粗细
            ax2.XAxis.LineWidth = 2;
            ax2.YAxis.LineWidth = 2;
            ax2.XAxis.FontSize = 20;
            ax2.YAxis.FontSize = 20;
            ax2.YAxis.FontWeight = 'bold';

            title_handle_2 = title('Cumulative Net Primary Productivity', 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
            ylabel('', 'FontName', 'Times New Roman','FontSize', 24, 'FontWeight', 'bold');
            
            if strcmp(season_file, 'Spring')
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 14, 'FontWeight', 'bold');
            else
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 18, 'FontWeight', 'bold');
            end
            h_legend.Box = 'off';
            
            % 获取坐标轴位置信息
            ax2_position2 = ax2.Position;

                                    % 调整 x 轴标签位置
            y_norm2 = 0.075/3 ;
            y_relative = y_norm2 * ax_position(4);

            abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            % 在原有刻度上添加额外的标记
            for i = 1:length(positions)
                x_position = positions(i);
                x_norm2 = (x_position - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));

                x_relative = x_norm2 * ax2_position2(3);

                % 绘制垂直于 x 轴的线条
                line([x_relative*(26/ax2_position2(3)), x_relative*(26/ax2_position2(3))], [abs_length + ylim_custom2(1),ylim_custom2(1)], ...
                    'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            
            % 调整 x 轴标签位置
            abs_length3 = y_relative3*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            ax2.XAxis.Label.Position(2) = ax2.XAxis.Label.Position(2) - abs_length3;

            % 循环添加标签
            for i = 1:length(positions)
                x_label = positions(i);
                x_norm_label = (x_label - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));
        
                x_relative_label = x_norm_label * ax2_position2(3);
        
                abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));
                
                % 添加标签
                text(x_relative_label*(26/ax2_position2(3)), ylim_custom2(1)-abs_length, labels{i}, 'Color', 'k', 'VerticalAlignment', 'top', ...
                    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
            end
            
            ax2.YAxis.Label.Position(1) = -1.5;
            
            % 调整垂直间距
            vertical_gap = 0.06;
            % 获取第一个子图的位置
            pos1 = get(h1, 'Position');
            
            % 获取第二个子图的位置，并调整其位置
            pos2 = get(h2, 'Position');
            pos2(2) = pos1(2) - pos1(4) + vertical_gap; % 调整第二个子图的 y 坐标位置
            
            % 更新第二个子图的位置
            set(h2, 'Position', pos2);

            box off     % 取消边框
            ax2_2 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax2_2,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax2_2.XAxis.LineWidth = 2;
            ax2_2.YAxis.LineWidth = 2;
            ax2_2.YAxis.FontWeight = 'bold';
            
            % 设置保存路径和文件名
            save_path = ['E:\Data\Drought_quantification\Quantify(pictures)\Quantify_0403\' segmented_name];
            
            % 完整的保存路径
            full_file_path = fullfile(save_path, output_file);
            
            % 设置输出分辨率为 300 DPI
            resolution = 500;
            
            % 保存图形为 JPEG 格式
            saveas(gcf, full_file_path, 'svg');
        elseif strcmp(data_kind, 'Ra')
            data_Ra = readmatrix(selected_file);
            data_Ra = data_Ra(:)';

            % 获取 period_points
            for n = 1:numel(selected_files_periods)
                selected_file_period = selected_files_periods{n};
                [~, period_kinds, ~] = fileparts(selected_file_period);
                period_words = split(period_kinds, '_');
                period_kind = period_words{4};
    
                if strcmp(period_kind, 'Ra') % 假设period_kind也应该是Ra
                    period_points = readmatrix(selected_file_period);
                    period_points = period_points(:)';
                    break; % 找到匹配的文件就退出循环
                end
            end

            Interpolation_values_drought = interp1(period, data_Ra, period_points);

            selected_file_path_parts = split(selected_file, filesep);
            fourth_file_name = selected_file_path_parts(5);

            if strcmp(fourth_file_name, 'Season_ALL_1_month')
                % 构造对应文件的路径
                corresponding_file = strrep(selected_file, 'Season_ALL_1_month', 'Season_ALL_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(10) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s.svg', data_kind, third_word, season_file);
            else
                corresponding_file = strrep(selected_file, 'Season_N&P_1_month', 'Season_N&P_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(11) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s(%s).svg', data_kind, third_word, season_file, tree_kind);
            end
            
            % 重新组合路径
            corresponding_file = fullfile(path_parts{:});

            data_Ra_no_drought = readmatrix(corresponding_file);
            data_Ra_no_drought = data_Ra_no_drought(:)';

            Interpolation_values_no_drought = interp1(period, data_Ra_no_drought, period_points);

            % 提取从 start_period 开始的时间和蒸腾变化量数据
            idx_start = find(period >= ceil(period_points(1)), 1);
            idx_end = find(period >= ceil(period_points(2)), 1);
            drought_period = [period_points(1), period(idx_start:idx_end-1),period_points(2)];
            drought_Ra_rate = [Interpolation_values_drought(1), data_Ra(idx_start:idx_end-1),Interpolation_values_drought(2)];
            
            % 干旱期间的累计蒸腾量
            cumulative_transpiration_droughting = cumtrapz(drought_period, drought_Ra_rate);
            
            new_drought_period = [drought_period, period(idx_end:end)];
            new_drought_Ra_rate = [drought_Ra_rate, data_Ra(idx_end:end)];
            
            cumulative_transpiration_drought = cumtrapz(new_drought_period, new_drought_Ra_rate);
            
            periods_before_drought = [period(1:idx_start - 1), period_points(1)];
            Ra_rate_before_drought = [data_Ra(1:idx_start - 1), Interpolation_values_drought(1)];
            cumulative_transpiration_before_drought = cumtrapz(periods_before_drought, Ra_rate_before_drought);
            
            new_period = [period(1:idx_start - 1), drought_period, period(idx_end:end)];
            new_Ra_rate = [data_Ra(1:idx_start - 1), drought_Ra_rate, data_Ra(idx_end:end)];
            
            % 假设没有发生干旱 假设没有发生干旱 假设没有发生干旱
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况

            integer_periods = ceil(period_points(1)):floor(period_points(2));
            
            % 计算整个时间范围内的累积蒸腾变化量
            cumulative_transpiration_total = cumtrapz(new_period, new_Ra_rate);
            
            lush_period = [period_points(1), integer_periods, period_points(2)];
            lush_Ra_rate = [Interpolation_values_no_drought(1), data_Ra_no_drought(integer_periods), Interpolation_values_no_drought(2)]; 
            
            % 计算没有干旱发生 的累计蒸腾量 针对干旱时期对应的时间
            cumulative_transpiration_lush = cumtrapz(lush_period, lush_Ra_rate);
            
            new_lush_period = [lush_period, period(idx_end:end)];
            new_lush_Ra_rate = [lush_Ra_rate, data_Ra_no_drought(idx_end:end)];
           
            % 计算没有干旱发生 的累计蒸腾量 针对干旱开始到25！！
            cumulative_transpiration_lush_after_drought = cumtrapz(new_lush_period, new_lush_Ra_rate);
            
            new_lush_period_total = [period(1:idx_start - 1), new_lush_period];
            new_lush_Ra_rate_total = [data_Ra_no_drought(1:idx_start - 1), new_lush_Ra_rate];
            
            cumulative_transpiration_lush_total = cumtrapz(new_lush_period_total, new_lush_Ra_rate_total);
            
            % 创建图形窗口并调整尺寸
            % figure('Units', 'normalized', 'Position', [0.6, 0.1, 0.3, 0.5]);
            figure('Units', 'pixels', 'Position', [1750,100,800,1000]);
            % 绘制蒸腾变化量随时间变化的曲线及其填充
            h1 = subplot(2, 1, 1);
            % 绘制干旱的情况
            p_star = plot(drought_period([1, end]), drought_Ra_rate([1, end]), 'p', 'MarkerSize', 3, 'LineWidth', 2 ,...
                'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255]);
            hold on;

%             p_original = plot(period, data_Ra, '-','Color', [255/255, 214/255, 0/255, 0.6],'LineWidth', 1.5);
%   
%             p_drought = plot(drought_period, drought_Ra_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.5); % 底部边界线
%             
%             p_no_drought = plot(period, data_Ra_no_drought, '-','Color', [255/255, 214/255, 0/255], 'LineWidth', 1.5);

            p_original = plot(period, data_Ra, '-k','LineWidth', 2);
            p_drought = plot(drought_period, drought_Ra_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5); % 底部边界线
            
            p_no_drought = plot(period, data_Ra_no_drought, '-','Color',[255/255, 214/255, 0/255], 'LineWidth', 2.5);

            xlim_custom = [0 26];        
            xticks_custom = 0:1:26;      
            ylim_custom = [-2 10];
            yticks_custom = -2:2:10;

            % 设置 x 和 y 轴范围
            xlim(xlim_custom);
            ylim(ylim_custom);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom);
            yticks(yticks_custom);

            set(gca, 'YTick', [0, 2, 4, 6, 8]);
            
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
            
            % 绘制新的时间点区间的蒸腾变化量曲线及其填充
            new_period_lush = [period(1:idx_start - 1), period_points(1), ceil(period_points(1)):floor(period_points(2)), period_points(2), period(idx_end:end)];
            new_Ra_rate_lush = [data_Ra(1:idx_start - 1), Interpolation_values_drought(1), interp1(period_points, Interpolation_values_drought, ceil(period_points(1)):floor(period_points(2))), Interpolation_values_drought(2), data_Ra(idx_end:end)];
            
            later_period = [period_points(2), period(idx_end:end)];
            later_Ra_rate = [Interpolation_values_drought(2), data_Ra(idx_end:end)];
            
            xverts_drought = [drought_period(1:end-1); drought_period(1:end-1); drought_period(2:end); drought_period(2:end)];
            yverts_drought = [zeros(1,length(drought_period)-1); drought_Ra_rate(1:end-1); drought_Ra_rate(2:end); zeros(1,length(drought_period)-1)];
            p_later = patch(xverts_drought, yverts_drought, [255/255, 214/255, 0/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            
            % 构造填充区域的顶点坐标
            xverts_fill = [drought_period, fliplr(lush_period)];
            yverts_fill = [drought_Ra_rate, fliplr(lush_Ra_rate)];
            
            % 绘制顶部边界线
            hold on;

            % 计算光带的上下边界
            offset = (0.08*9)/16; % 定义偏移量，可以调整这个值来改变光带的宽度
            upper_bound = data_Ra_no_drought + offset;
            lower_bound = data_Ra_no_drought - offset;

%             % 第一层光带，透明度较高
%             fill([period, fliplr(period)], [upper_bound+(0.08*9)/16, fliplr(lower_bound-(0.08*9)/16)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
%             % 第二层光带，透明度较低
%             fill([period, fliplr(period)], [upper_bound, fliplr(lower_bound)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

            % 绘制不干旱的情况
            plot(period, data_Ra_no_drought, '-','Color',[255/255, 214/255, 0/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');

            % 使用亮眼的颜色进行填充，并设置透明度
            % highlight_color = [184/255, 160/255, 52/255];%枯萎！！
            highlight_color = [211/255, 63/255, 43/255];
            p_fill = patch(xverts_fill, yverts_fill, highlight_color, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 无边界颜色

            plot(drought_period, drought_Ra_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off'); % 底部边界线

            plot(drought_period([1, end]), drought_Ra_rate([1, end]), 'p', ...
                'MarkerSize', 3, 'LineWidth', 2, 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255], 'HandleVisibility', 'off');
            
            % 设置标题和 y 轴标签，根据 third_word 和 drought_speed
            if isempty(tree_kind)
                if strcmp(third_word, 'FD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'Ra')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = '';
                    end
                end
            else
                if strcmp(third_word, 'FD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'Ra')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'Ra')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    end
                end
            end

            % 获取当前坐标轴对象
            ax = gca;
            set(gca,  'XTickLabel', []);
            
            % 获取坐标轴位置信息
            ax_position = ax.Position;
            
            ax.YAxis.Label.Position(1) = -1.5;
            
            % 设置 x 轴和 y 轴的线条粗细
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            ax.YAxis.FontWeight = 'bold';

                                                % 获取当前子图的当前位置
            pos = get(gca, 'Position');
            
            % 增加高度并保持上端位置不变
            new_height = pos(4) * 1.3;  % 增加高度
            pos(2) = pos(2) - (new_height - pos(4));  % 调整下端位置
            pos(4) = new_height;  % 设置新的高度
            
            % 应用新的位置
            set(gca, 'Position', pos);

            title_handle_1 = title(title_str, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            ylabel(ylabel_str, 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');

            if strcmp(season_file, 'Spring')
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Respiration', 'Respiration Anomaly'}, 'Location', 'northwest','FontSize', 14, 'FontWeight', 'bold');
            else
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Respiration', 'Respiration Anomaly'}, 'Location', 'northeast','FontSize', 18, 'FontWeight', 'bold');
            end
            legend_handle.Box = 'off';

            y_norm2 = 0.06/3;
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
            y_norm3 = 0.15/3;
            y_relative3 = y_norm3 * ax_position(4);

            box off     % 取消边框
            ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax1.XAxis.LineWidth = 2;
            ax1.YAxis.LineWidth = 2;
            ax1.YAxis.FontWeight = 'bold';

            

            % 第二个子图：绘制新的时间点区间的蒸腾变化量曲线和填充色
            h2 = subplot(2, 1, 2);
            plot(periods_before_drought, cumulative_transpiration_before_drought, '-k','LineWidth', 2.5);
            hold on
            
            offset = cumulative_transpiration_before_drought(end); % 你想加上的值
            
            % 修改后的数据
            cumulative_transpiration_drought_offset = cumulative_transpiration_droughting + offset;
            cumulative_transpiration_lush_offset = cumulative_transpiration_lush + offset;
            
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5);
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [255/255, 214/255, 0/255],'LineWidth',2.5);

            % 根据 fourth_word 设置不同的 x 和 y 轴范围和刻度间隔
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 75];    
                            yticks_custom2 = 0:15:75;
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 75];    
                            yticks_custom2 = 0:15:75;
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        xlim_custom2 = [0 26];        
                        xticks_custom2 = 0:1:26;      
                        ylim_custom2 = [0 75];    
                        yticks_custom2 = 0:15:75;
                end
            end

            % 设置 x 和 y 轴范围
            xlim(xlim_custom2);
            ylim(ylim_custom2);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom2);
            yticks(yticks_custom2);
            
            % 创建阴影区域
            x_fill_shadow_1 = [period_points(1), serious_period, serious_period, period_points(1)];
            y_fill_shadow_1 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_1, y_fill_shadow_1, [128/255, 128/255, 128/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_2 = [serious_period, period_points(2), period_points(2), serious_period];
            y_fill_shadow_2 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_2, y_fill_shadow_2, [128/255, 128/255, 128/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_3 = [period_points(2), 26, 26, period_points(2)];
            y_fill_shadow_3 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_3, y_fill_shadow_3, [128/255, 128/255, 128/255], 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            hold on
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [255/255, 214/255, 0/255],'LineWidth',2.5, 'HandleVisibility', 'off');
            
            % 添加箭头、短线
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file

                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 214/255, 0/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 214/255, 0/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.8, '^r', ...
                            'MarkerFaceColor', [255/255, 214/255, 0/255], 'MarkerEdgeColor', [255/255, 214/255, 0/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+1.6, 'vr', ...
                            'MarkerFaceColor', [255/255, 214/255, 0/255], 'MarkerEdgeColor', [255/255, 214/255, 0/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [255/255, 214/255, 0/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-1.6, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.8, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);
                    end

                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        
                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 214/255, 0/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 214/255, 0/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-1, '^r', ...
                            'MarkerFaceColor', [255/255, 214/255, 0/255], 'MarkerEdgeColor', [255/255, 214/255, 0/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+1.5, 'vr', ...
                            'MarkerFaceColor', [255/255, 214/255, 0/255], 'MarkerEdgeColor', [255/255, 214/255, 0/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [255/255, 214/255, 0/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-1.5, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+1, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);

                    end
                end
            else
                switch season_file

                    case 'Autumn'
                    % 对于不干旱的情况
                    % 绘制横线
                    line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 214/255, 0/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 214/255, 0/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(lush_period(1), cumulative_transpiration_lush_offset(end)-1, '^r', ...
                        'MarkerFaceColor', [255/255, 214/255, 0/255], 'MarkerEdgeColor', [255/255, 214/255, 0/255], 'MarkerSize', 5);  % 上箭头
                    plot(lush_period(1), cumulative_transpiration_lush_offset(1)+1.5, 'vr', ...
                        'MarkerFaceColor', [255/255, 214/255, 0/255], 'MarkerEdgeColor', [255/255, 214/255, 0/255],'MarkerSize', 5);  % 下箭头
                    
                    % 计算线段长度
                    length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                        ['{\bf ', num2str(round(length_of_line, 2)), '} \bf{}'], ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [255/255, 214/255, 0/255]);
                    
                    % 对于干旱情况
                    % 绘制横线
                    line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(drought_period(end), cumulative_transpiration_drought_offset(end)-1.5, '^r', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                    plot(drought_period(end), cumulative_transpiration_drought_offset(1)+1, 'vr', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                    
                    
                    % 计算线段长度
                    length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                        ['{\bf}', num2str(round(length_of_line2, 2))], ...
                        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);

                end
            end
            
            scatter(drought_period(1), cumulative_transpiration_drought_offset(1), ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(lush_period(end), cumulative_transpiration_lush_offset(end), ...
                'MarkerFaceColor', [255/255, 214/255, 0/255], 'MarkerEdgeColor', [255/255, 214/255, 0/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            % grid on;
            % 黑色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.6, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.6, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.7, ...
                            ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                            'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                end
            end

            % 红色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.04, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.02, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                end
            end

            % 绿色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(lush_period(end), cumulative_transpiration_lush_offset(end)+3.5, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [255/255, 214/255, 0/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                         text(lush_period(end), cumulative_transpiration_lush_offset(end)+3, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [255/255, 214/255, 0/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                     text(lush_period(end), cumulative_transpiration_lush_offset(end)+3, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [255/255, 214/255, 0/255]);
                end
            end

            
            % 获取当前坐标轴对象
            ax2 = gca;
            set(gca,  'XTickLabel', []);

            % 设置 x 轴和 y 轴的线条粗细
            ax2.XAxis.LineWidth = 2;
            ax2.YAxis.LineWidth = 2;
            ax2.XAxis.FontSize = 20;
            ax2.YAxis.FontSize = 20;
            ax2.YAxis.FontWeight = 'bold';

            title_handle_2 = title('Cumulative Autotrophic Respiration', 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
            ylabel('', 'FontName', 'Times New Roman','FontSize', 24, 'FontWeight', 'bold');
            
            if strcmp(season_file, 'Spring')
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 14, 'FontWeight', 'bold');
            else
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 18, 'FontWeight', 'bold');
            end
            h_legend.Box = 'off';
            
            % 获取坐标轴位置信息
            ax2_position2 = ax2.Position;

                                    % 调整 x 轴标签位置
            y_norm2 = 0.075/3 ;
            y_relative = y_norm2 * ax_position(4);

            abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            % 在原有刻度上添加额外的标记
            for i = 1:length(positions)
                x_position = positions(i);
                x_norm2 = (x_position - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));

                x_relative = x_norm2 * ax2_position2(3);

                % 绘制垂直于 x 轴的线条
                line([x_relative*(26/ax2_position2(3)), x_relative*(26/ax2_position2(3))], [abs_length + ylim_custom2(1),ylim_custom2(1)], ...
                    'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            
            % 调整 x 轴标签位置
            abs_length3 = y_relative3*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            ax2.XAxis.Label.Position(2) = ax2.XAxis.Label.Position(2) - abs_length3;

            % 循环添加标签
            for i = 1:length(positions)
                x_label = positions(i);
                x_norm_label = (x_label - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));
        
                x_relative_label = x_norm_label * ax2_position2(3);
        
                abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));
                
                % 添加标签
                text(x_relative_label*(26/ax2_position2(3)), ylim_custom2(1)-abs_length, labels{i}, 'Color', 'k', 'VerticalAlignment', 'top', ...
                    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
            end
            
            ax2.YAxis.Label.Position(1) = -1.5;
            
            % 调整垂直间距
            vertical_gap = 0.06;
            % 获取第一个子图的位置
            pos1 = get(h1, 'Position');
            
            % 获取第二个子图的位置，并调整其位置
            pos2 = get(h2, 'Position');
            pos2(2) = pos1(2) - pos1(4) + vertical_gap; % 调整第二个子图的 y 坐标位置
            
            % 更新第二个子图的位置
            set(h2, 'Position', pos2);

            box off     % 取消边框
            ax2_2 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax2_2,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax2_2.XAxis.LineWidth = 2;
            ax2_2.YAxis.LineWidth = 2;
            ax2_2.YAxis.FontWeight = 'bold';
            
            % 设置保存路径和文件名
            save_path = ['E:\Data\Drought_quantification\Quantify(pictures)\Quantify_0403\' segmented_name];
            
            % 完整的保存路径
            full_file_path = fullfile(save_path, output_file);
            
            % 设置输出分辨率为 300 DPI
            resolution = 500;
            
            % 保存图形为 JPEG 格式
            saveas(gcf, full_file_path, 'svg');

        elseif strcmp(data_kind, 'Reco')
            data_Reco = readmatrix(selected_file);
            data_Reco = data_Reco(:)';

            % 获取 period_points
            for n = 1:numel(selected_files_periods)
                selected_file_period = selected_files_periods{n};
                [~, period_kinds, ~] = fileparts(selected_file_period);
                period_words = split(period_kinds, '_');
                period_kind = period_words{4};
    
                if strcmp(period_kind, 'Reco') % 假设period_kind也应该是Reco
                    period_points = readmatrix(selected_file_period);
                    period_points = period_points(:)';
                    break; % 找到匹配的文件就退出循环
                end
            end

            Interpolation_values_drought = interp1(period, data_Reco, period_points);

            selected_file_path_parts = split(selected_file, filesep);
            fourth_file_name = selected_file_path_parts(5);

            if strcmp(fourth_file_name, 'Season_ALL_1_month')
                % 构造对应文件的路径
                corresponding_file = strrep(selected_file, 'Season_ALL_1_month', 'Season_ALL_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(10) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s.svg', data_kind, third_word, season_file);
            else
                corresponding_file = strrep(selected_file, 'Season_N&P_1_month', 'Season_N&P_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(11) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s(%s).svg', data_kind, third_word, season_file, tree_kind);
            end
            
            % 重新组合路径
            corresponding_file = fullfile(path_parts{:});

            data_Reco_no_drought = readmatrix(corresponding_file);
            data_Reco_no_drought = data_Reco_no_drought(:)';

            Interpolation_values_no_drought = interp1(period, data_Reco_no_drought, period_points);

            % 提取从 start_period 开始的时间和蒸腾变化量数据
            idx_start = find(period >= ceil(period_points(1)), 1);
            idx_end = find(period >= ceil(period_points(2)), 1);
            drought_period = [period_points(1), period(idx_start:idx_end-1),period_points(2)];
            drought_Reco_rate = [Interpolation_values_drought(1), data_Reco(idx_start:idx_end-1),Interpolation_values_drought(2)];
            
            % 干旱期间的累计蒸腾量
            cumulative_transpiration_droughting = cumtrapz(drought_period, drought_Reco_rate);
            
            new_drought_period = [drought_period, period(idx_end:end)];
            new_drought_Reco_rate = [drought_Reco_rate, data_Reco(idx_end:end)];
            
            cumulative_transpiration_drought = cumtrapz(new_drought_period, new_drought_Reco_rate);
            
            periods_before_drought = [period(1:idx_start - 1), period_points(1)];
            Reco_rate_before_drought = [data_Reco(1:idx_start - 1), Interpolation_values_drought(1)];
            cumulative_transpiration_before_drought = cumtrapz(periods_before_drought, Reco_rate_before_drought);
            
            new_period = [period(1:idx_start - 1), drought_period, period(idx_end:end)];
            new_Reco_rate = [data_Reco(1:idx_start - 1), drought_Reco_rate, data_Reco(idx_end:end)];
            
            % 假设没有发生干旱 假设没有发生干旱 假设没有发生干旱
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况

            integer_periods = ceil(period_points(1)):floor(period_points(2));
            
            % 计算整个时间范围内的累积蒸腾变化量
            cumulative_transpiration_total = cumtrapz(new_period, new_Reco_rate);
            
            lush_period = [period_points(1), integer_periods, period_points(2)];
            lush_Reco_rate = [Interpolation_values_no_drought(1), data_Reco_no_drought(integer_periods), Interpolation_values_no_drought(2)]; 
            
            % 计算没有干旱发生 的累计蒸腾量 针对干旱时期对应的时间
            cumulative_transpiration_lush = cumtrapz(lush_period, lush_Reco_rate);
            
            new_lush_period = [lush_period, period(idx_end:end)];
            new_lush_Reco_rate = [lush_Reco_rate, data_Reco_no_drought(idx_end:end)];
           
            % 计算没有干旱发生 的累计蒸腾量 针对干旱开始到25！！
            cumulative_transpiration_lush_after_drought = cumtrapz(new_lush_period, new_lush_Reco_rate);
            
            new_lush_period_total = [period(1:idx_start - 1), new_lush_period];
            new_lush_Reco_rate_total = [data_Reco_no_drought(1:idx_start - 1), new_lush_Reco_rate];
            
            cumulative_transpiration_lush_total = cumtrapz(new_lush_period_total, new_lush_Reco_rate_total);
            
            % 创建图形窗口并调整尺寸
            % figure('Units', 'normalized', 'Position', [0.6, 0.1, 0.3, 0.5]);
            figure('Units', 'pixels', 'Position', [1750,100,800,1000]);
            % 绘制蒸腾变化量随时间变化的曲线及其填充
            h1 = subplot(2, 1, 1);
            % 绘制干旱的情况
            p_star = plot(drought_period([1, end]), drought_Reco_rate([1, end]), 'p', 'MarkerSize', 3, 'LineWidth', 2 ,...
                'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255]);
            hold on;
%             p_original = plot(period, data_Reco, '-','Color', [255/255, 120/255, 13/255, 0.8],'LineWidth', 1.5);
%   
%             p_drought = plot(drought_period, drought_Reco_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.5); % 底部边界线
%             
%             p_no_drought = plot(period, data_Reco_no_drought, '-','Color', [255/255, 120/255, 13/255], 'LineWidth', 1.5);

            p_original = plot(period, data_Reco, '-k','LineWidth', 2);
            p_drought = plot(drought_period, drought_Reco_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5); % 底部边界线
            
            p_no_drought = plot(period, data_Reco_no_drought, '-','Color',[255/255, 120/255, 13/255], 'LineWidth', 2.5);

            xlim_custom = [0 26];        
            xticks_custom = 0:1:26;      
            ylim_custom = [-2 12];    
            yticks_custom = -2:3:12;

            % 设置 x 和 y 轴范围
            xlim(xlim_custom);
            ylim(ylim_custom);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom);
            yticks(yticks_custom);

            set(gca, 'YTick', [0, 2, 4, 6, 8, 10]);
            
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
            
            % 绘制新的时间点区间的蒸腾变化量曲线及其填充
            new_period_lush = [period(1:idx_start - 1), period_points(1), ceil(period_points(1)):floor(period_points(2)), period_points(2), period(idx_end:end)];
            new_Reco_rate_lush = [data_Reco(1:idx_start - 1), Interpolation_values_drought(1), interp1(period_points, Interpolation_values_drought, ceil(period_points(1)):floor(period_points(2))), Interpolation_values_drought(2), data_Reco(idx_end:end)];
            
            later_period = [period_points(2), period(idx_end:end)];
            later_Reco_rate = [Interpolation_values_drought(2), data_Reco(idx_end:end)];
            
            xverts_drought = [drought_period(1:end-1); drought_period(1:end-1); drought_period(2:end); drought_period(2:end)];
            yverts_drought = [zeros(1,length(drought_period)-1); drought_Reco_rate(1:end-1); drought_Reco_rate(2:end); zeros(1,length(drought_period)-1)];
            p_later = patch(xverts_drought, yverts_drought, [255/255, 120/255, 13/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            
            % 构造填充区域的顶点坐标
            xverts_fill = [drought_period, fliplr(lush_period)];
            yverts_fill = [drought_Reco_rate, fliplr(lush_Reco_rate)];
            
            % 绘制顶部边界线
            hold on;

            % 计算光带的上下边界
            offset = (0.08*15)/16; % 定义偏移量，可以调整这个值来改变光带的宽度
            upper_bound = data_Reco_no_drought + offset;
            lower_bound = data_Reco_no_drought - offset;

%             % 第一层光带，透明度较高
%             fill([period, fliplr(period)], [upper_bound+(0.08*15)/16, fliplr(lower_bound-(0.08*15)/16)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
%             % 第二层光带，透明度较低
%             fill([period, fliplr(period)], [upper_bound, fliplr(lower_bound)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

            % 绘制不干旱的情况
            plot(period, data_Reco_no_drought, '-','Color',[255/255, 120/255, 13/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');

            % 使用亮眼的颜色进行填充，并设置透明度
            % highlight_color = [184/255, 160/255, 52/255];%枯萎！！
            highlight_color = [211/255, 63/255, 43/255];
            p_fill = patch(xverts_fill, yverts_fill, highlight_color, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 无边界颜色

            plot(drought_period, drought_Reco_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off'); % 底部边界线

            plot(drought_period([1, end]), drought_Reco_rate([1, end]), 'p', ...
                'MarkerSize', 3, 'LineWidth', 2, 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255], 'HandleVisibility', 'off');
            
            % 设置标题和 y 轴标签，根据 third_word 和 drought_speed
            if isempty(tree_kind)
                if strcmp(third_word, 'FD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'Reco')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = '';
                    end
                end
            else
                if strcmp(third_word, 'FD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'Reco')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'Reco')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    end
                end
            end

            % 获取当前坐标轴对象
            ax = gca;
            set(gca,  'XTickLabel', []);
            
            % 获取坐标轴位置信息
            ax_position = ax.Position;
            
            ax.YAxis.Label.Position(1) = -1.5;
            
            % 设置 x 轴和 y 轴的线条粗细
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            ax.YAxis.FontWeight = 'bold';

                                                % 获取当前子图的当前位置
            pos = get(gca, 'Position');
            
            % 增加高度并保持上端位置不变
            new_height = pos(4) * 1.3;  % 增加高度
            pos(2) = pos(2) - (new_height - pos(4));  % 调整下端位置
            pos(4) = new_height;  % 设置新的高度
            
            % 应用新的位置
            set(gca, 'Position', pos);

            title_handle_1 = title(title_str, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            ylabel(ylabel_str, 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');

            if strcmp(season_file, 'Spring')
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Respiration', 'Respiration Anomaly'}, 'Location', 'northwest','FontSize', 14, 'FontWeight', 'bold');
            else
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Respiration', 'Respiration Anomaly'}, 'Location', 'northeast','FontSize', 18, 'FontWeight', 'bold');
            end
            legend_handle.Box = 'off';

            y_norm2 = 0.06/3;
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
            y_norm3 = 0.15/3;
            y_relative3 = y_norm3 * ax_position(4);

            box off     % 取消边框
            ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax1.XAxis.LineWidth = 2;
            ax1.YAxis.LineWidth = 2;
            ax1.YAxis.FontWeight = 'bold';

            

            % 第二个子图：绘制新的时间点区间的蒸腾变化量曲线和填充色
            h2 = subplot(2, 1, 2);
            plot(periods_before_drought, cumulative_transpiration_before_drought, '-k','LineWidth', 2.5);
            hold on
            
            offset = cumulative_transpiration_before_drought(end); % 你想加上的值
            
            % 修改后的数据
            cumulative_transpiration_drought_offset = cumulative_transpiration_droughting + offset;
            cumulative_transpiration_lush_offset = cumulative_transpiration_lush + offset;
            
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5);
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [255/255, 120/255, 13/255],'LineWidth',2.5);

            % 根据 fourth_word 设置不同的 x 和 y 轴范围和刻度间隔
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 120];    
                            yticks_custom2 = 0:30:120;
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 120];    
                            yticks_custom2 = 0:30:120;
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        xlim_custom2 = [0 26];        
                        xticks_custom2 = 0:1:26;      
                        ylim_custom2 = [0 120];    
                        yticks_custom2 = 0:30:120;
                end
            end

            % 设置 x 和 y 轴范围
            xlim(xlim_custom2);
            ylim(ylim_custom2);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom2);
            yticks(yticks_custom2);
            
            % 创建阴影区域
            x_fill_shadow_1 = [period_points(1), serious_period, serious_period, period_points(1)];
            y_fill_shadow_1 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_1, y_fill_shadow_1, [128/255, 128/255, 128/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_2 = [serious_period, period_points(2), period_points(2), serious_period];
            y_fill_shadow_2 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_2, y_fill_shadow_2, [128/255, 128/255, 128/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_3 = [period_points(2), 26, 26, period_points(2)];
            y_fill_shadow_3 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_3, y_fill_shadow_3, [128/255, 128/255, 128/255], 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            hold on
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [255/255, 120/255, 13/255],'LineWidth',2.5, 'HandleVisibility', 'off');
            
            % 添加箭头、短线
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file

                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 120/255, 13/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 120/255, 13/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-1.5, '^r', ...
                            'MarkerFaceColor', [255/255, 120/255, 13/255], 'MarkerEdgeColor', [255/255, 120/255, 13/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+3, 'vr', ...
                            'MarkerFaceColor', [255/255, 120/255, 13/255], 'MarkerEdgeColor', [255/255, 120/255, 13/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [255/255, 120/255, 13/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-3, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+1.5, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);
                    end

                elseif strcmp(tree_kind,'Planted')
                    switch season_file

                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 120/255, 13/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 120/255, 13/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-1.5, '^r', ...
                            'MarkerFaceColor', [255/255, 120/255, 13/255], 'MarkerEdgeColor', [255/255, 120/255, 13/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+3, 'vr', ...
                            'MarkerFaceColor', [255/255, 120/255, 13/255], 'MarkerEdgeColor', [255/255, 120/255, 13/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [255/255, 120/255, 13/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-3, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+1.5, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    % 对于不干旱的情况
                    % 绘制横线
                    line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 120/255, 13/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [255/255, 120/255, 13/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(lush_period(1), cumulative_transpiration_lush_offset(end)-1.5, '^r', ...
                        'MarkerFaceColor', [255/255, 120/255, 13/255], 'MarkerEdgeColor', [255/255, 120/255, 13/255], 'MarkerSize', 5);  % 上箭头
                    plot(lush_period(1), cumulative_transpiration_lush_offset(1)+3, 'vr', ...
                        'MarkerFaceColor', [255/255, 120/255, 13/255], 'MarkerEdgeColor', [255/255, 120/255, 13/255],'MarkerSize', 5);  % 下箭头
                    
                    % 计算线段长度
                    length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                        ['{\bf ', num2str(round(length_of_line, 2)), '} \bf{}'], ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [255/255, 120/255, 13/255]);
                    
                    % 对于干旱情况
                    % 绘制横线
                    line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(drought_period(end), cumulative_transpiration_drought_offset(end)-3, '^r', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                    plot(drought_period(end), cumulative_transpiration_drought_offset(1)+1, 'vr', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                    
                    
                    % 计算线段长度
                    length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                    
                    % 添加长度标注并在文本中添加短线
                    text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                        ['{\bf}', num2str(round(length_of_line2, 2))], ...
                        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                end
            end
            
            scatter(drought_period(1), cumulative_transpiration_drought_offset(1), ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(lush_period(end), cumulative_transpiration_lush_offset(end), ...
                'MarkerFaceColor', [255/255, 120/255, 13/255], 'MarkerEdgeColor', [255/255, 120/255, 13/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            % grid on;
            % 黑色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+1, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+1, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.05, ...
                            ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{—}'], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                            'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                end
            end

            % 红色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.04, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.02, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                end
            end

            % 绿色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(lush_period(end), cumulative_transpiration_lush_offset(end)+6, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [255/255, 120/255, 13/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                         text(lush_period(end), cumulative_transpiration_lush_offset(end)+6, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [255/255, 120/255, 13/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                     text(lush_period(end), cumulative_transpiration_lush_offset(end)+5, ...
                        ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [255/255, 120/255, 13/255]);
                end
            end

            
            % 获取当前坐标轴对象
            ax2 = gca;
            set(gca,  'XTickLabel', []);

            % 设置 x 轴和 y 轴的线条粗细
            ax2.XAxis.LineWidth = 2;
            ax2.YAxis.LineWidth = 2;
            ax2.XAxis.FontSize = 20;
            ax2.YAxis.FontSize = 20;
            ax2.YAxis.FontWeight = 'bold';

            title_handle_2 = title('Cumulative Ecosystem Respiration', 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
            ylabel('', 'FontName', 'Times New Roman','FontSize', 24, 'FontWeight', 'bold');
            
            if strcmp(season_file, 'Spring')
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 14, 'FontWeight', 'bold');
            else
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 18, 'FontWeight', 'bold');
            end
            h_legend.Box = 'off';
            
            % 获取坐标轴位置信息
            ax2_position2 = ax2.Position;

                                    % 调整 x 轴标签位置
            y_norm2 = 0.075/3 ;
            y_relative = y_norm2 * ax_position(4);

            abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            % 在原有刻度上添加额外的标记
            for i = 1:length(positions)
                x_position = positions(i);
                x_norm2 = (x_position - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));

                x_relative = x_norm2 * ax2_position2(3);

                % 绘制垂直于 x 轴的线条
                line([x_relative*(26/ax2_position2(3)), x_relative*(26/ax2_position2(3))], [abs_length + ylim_custom2(1),ylim_custom2(1)], ...
                    'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            
            % 调整 x 轴标签位置
            abs_length3 = y_relative3*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            ax2.XAxis.Label.Position(2) = ax2.XAxis.Label.Position(2) - abs_length3;

            % 循环添加标签
            for i = 1:length(positions)
                x_label = positions(i);
                x_norm_label = (x_label - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));
        
                x_relative_label = x_norm_label * ax2_position2(3);
        
                abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));
                
                % 添加标签
                text(x_relative_label*(26/ax2_position2(3)), ylim_custom2(1)-abs_length, labels{i}, 'Color', 'k', 'VerticalAlignment', 'top', ...
                    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
            end
            
            ax2.YAxis.Label.Position(1) = -1.5;
            
            % 调整垂直间距
            vertical_gap = 0.06;
            % 获取第一个子图的位置
            pos1 = get(h1, 'Position');
            
            % 获取第二个子图的位置，并调整其位置
            pos2 = get(h2, 'Position');
            pos2(2) = pos1(2) - pos1(4) + vertical_gap; % 调整第二个子图的 y 坐标位置
            
            % 更新第二个子图的位置
            set(h2, 'Position', pos2);

            box off     % 取消边框
            ax2_2 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax2_2,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax2_2.XAxis.LineWidth = 2;
            ax2_2.YAxis.LineWidth = 2;
            ax2_2.YAxis.FontWeight = 'bold';
            
            % 设置保存路径和文件名
            save_path = ['E:\Data\Drought_quantification\Quantify(pictures)\Quantify_0403\' segmented_name];
            
            % 完整的保存路径
            full_file_path = fullfile(save_path, output_file);
            
            % 设置输出分辨率为 300 DPI
            resolution = 500;
            
            % 保存图形为 JPEG 格式
            saveas(gcf, full_file_path, 'svg');

        elseif strcmp(data_kind,'T')
            data_T = readmatrix(selected_file);
            data_T = data_T(:)';

            % 获取 period_points
            for n = 1:numel(selected_files_periods)
                selected_file_period = selected_files_periods{n};
                [~, period_kinds, ~] = fileparts(selected_file_period);
                period_words = split(period_kinds, '_');
                period_kind = period_words{4};
    
                if strcmp(period_kind, 'T') % 假设period_kind也应该是T
                    period_points = readmatrix(selected_file_period);
                    period_points = period_points(:)';
                    break; % 找到匹配的文件就退出循环
                end
            end

            Interpolation_values_drought = interp1(period, data_T, period_points);

            selected_file_path_parts = split(selected_file, filesep);
            fourth_file_name = selected_file_path_parts(5);

            if strcmp(fourth_file_name, 'Season_ALL_1_month')
                % 构造对应文件的路径
                corresponding_file = strrep(selected_file, 'Season_ALL_1_month', 'Season_ALL_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(10) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s.svg', data_kind, third_word, season_file);
            else
                corresponding_file = strrep(selected_file, 'Season_N&P_1_month', 'Season_N&P_no_drought_1_month');

                % 分割路径为目录和文件名部分
                path_parts = split(corresponding_file, filesep);
                
                % 去掉需要删除的目录
                path_parts(11) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s(%s).svg', data_kind, third_word, season_file, tree_kind);
            end
            
            % 重新组合路径
            corresponding_file = fullfile(path_parts{:});

            data_T_no_drought = readmatrix(corresponding_file);
            data_T_no_drought = data_T_no_drought(:)';

            Interpolation_values_no_drought = interp1(period, data_T_no_drought, period_points);

            % 提取从 start_period 开始的时间和蒸腾变化量数据
            idx_start = find(period >= ceil(period_points(1)), 1);
            idx_end = find(period >= ceil(period_points(2)), 1);
            drought_period = [period_points(1), period(idx_start:idx_end-1),period_points(2)];
            drought_T_rate = [Interpolation_values_drought(1), data_T(idx_start:idx_end-1),Interpolation_values_drought(2)];
            
            % 干旱期间的累计蒸腾量
            cumulative_transpiration_droughting = cumtrapz(drought_period, drought_T_rate);
            
            new_drought_period = [drought_period, period(idx_end:end)];
            new_drought_T_rate = [drought_T_rate, data_T(idx_end:end)];
            
            cumulative_transpiration_drought = cumtrapz(new_drought_period, new_drought_T_rate);
            
            periods_before_drought = [period(1:idx_start - 1), period_points(1)];
            T_rate_before_drought = [data_T(1:idx_start - 1), Interpolation_values_drought(1)];
            cumulative_transpiration_before_drought = cumtrapz(periods_before_drought, T_rate_before_drought);
            
            new_period = [period(1:idx_start - 1), drought_period, period(idx_end:end)];
            new_T_rate = [data_T(1:idx_start - 1), drought_T_rate, data_T(idx_end:end)];
            
            % 假设没有发生干旱 假设没有发生干旱 假设没有发生干旱
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况
            % 计算如果没有发生干旱的情况

            integer_periods = ceil(period_points(1)):floor(period_points(2));
            
            % 计算整个时间范围内的累积蒸腾变化量
            cumulative_transpiration_total = cumtrapz(new_period, new_T_rate);
            
            lush_period = [period_points(1), integer_periods, period_points(2)];
            lush_T_rate = [Interpolation_values_no_drought(1), data_T_no_drought(integer_periods), Interpolation_values_no_drought(2)]; 
            
            % 计算没有干旱发生 的累计蒸腾量 针对干旱时期对应的时间
            cumulative_transpiration_lush = cumtrapz(lush_period, lush_T_rate);
            
            new_lush_period = [lush_period, period(idx_end:end)];
            new_lush_T_rate = [lush_T_rate, data_T_no_drought(idx_end:end)];
           
            % 计算没有干旱发生 的累计蒸腾量 针对干旱开始到25！！
            cumulative_transpiration_lush_after_drought = cumtrapz(new_lush_period, new_lush_T_rate);
            
            new_lush_period_total = [period(1:idx_start - 1), new_lush_period];
            new_lush_T_rate_total = [data_T_no_drought(1:idx_start - 1), new_lush_T_rate];
            
            cumulative_transpiration_lush_total = cumtrapz(new_lush_period_total, new_lush_T_rate_total);
            
            % 创建图形窗口并调整尺寸
            % figure('Units', 'normalized', 'Position', [0.6, 0.1, 0.3, 0.5]);
            figure('Units', 'pixels', 'Position', [1750,100,800,1000]);
            % 绘制蒸腾变化量随时间变化的曲线及其填充
            h1 = subplot(2, 1, 1);
            % 绘制干旱的情况
            p_star = plot(drought_period([1, end]), drought_T_rate([1, end]), 'p', 'MarkerSize', 3, 'LineWidth', 2 ,...
                'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255]);
            hold on;

%             p_original = plot(period, data_T, '-','Color', [12/255, 114/255, 31/255, 0.8],'LineWidth', 1.5);
%   
%             p_drought = plot(drought_period, drought_T_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.5); % 底部边界线
%             
%             p_no_drought = plot(period, data_T_no_drought, '-','Color', [12/255, 114/255, 31/255], 'LineWidth', 1.5);

            p_original = plot(period, data_T, '-k','LineWidth', 2);
            p_drought = plot(drought_period, drought_T_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5); % 底部边界线
            
            p_no_drought = plot(period, data_T_no_drought, '-','Color',[12/255, 114/255, 31/255], 'LineWidth', 2.5);

            xlim_custom = [0 26];        
            xticks_custom = 0:1:26;      
            ylim_custom = [-0.8 4.8]; 
%             yticks_custom = -1:2:5;

            % 设置 x 和 y 轴范围
            xlim(xlim_custom);
            ylim(ylim_custom);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom);
            set(gca, 'YTick', [0, 0.8, 1.6, 2.4, 3.2, 4.0]); % 设置 y 轴的刻度为 -1, 0, 2, 4, 6);

            
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
            
            % 绘制新的时间点区间的蒸腾变化量曲线及其填充
            new_period_lush = [period(1:idx_start - 1), period_points(1), ceil(period_points(1)):floor(period_points(2)), period_points(2), period(idx_end:end)];
            new_T_rate_lush = [data_T(1:idx_start - 1), Interpolation_values_drought(1), interp1(period_points, Interpolation_values_drought, ceil(period_points(1)):floor(period_points(2))), Interpolation_values_drought(2), data_T(idx_end:end)];
            
            later_period = [period_points(2), period(idx_end:end)];
            later_T_rate = [Interpolation_values_drought(2), data_T(idx_end:end)];
            
            xverts_drought = [drought_period(1:end-1); drought_period(1:end-1); drought_period(2:end); drought_period(2:end)];
            yverts_drought = [zeros(1,length(drought_period)-1); drought_T_rate(1:end-1); drought_T_rate(2:end); zeros(1,length(drought_period)-1)];
            p_later = patch(xverts_drought, yverts_drought, [12/255, 114/255, 31/255], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
            
            % 构造填充区域的顶点坐标
            xverts_fill = [drought_period, fliplr(lush_period)];
            yverts_fill = [drought_T_rate, fliplr(lush_T_rate)];
            
            % 绘制顶部边界线
            hold on;

            % 计算光带的上下边界
            offset = (0.08*7)/16; % 定义偏移量，可以调整这个值来改变光带的宽度
            upper_bound = data_T_no_drought + offset;
            lower_bound = data_T_no_drought - offset;

            % 第一层光带，透明度较高
%             fill([period, fliplr(period)], [upper_bound+(0.08*7)/16, fliplr(lower_bound-(0.08*7)/16)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
%             % 第二层光带，透明度较低
%             fill([period, fliplr(period)], [upper_bound, fliplr(lower_bound)], [151/255, 253/255, 82/255], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

            % 绘制不干旱的情况
            plot(period, data_T_no_drought, '-','Color',[12/255, 114/255, 31/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');

            % 使用亮眼的颜色进行填充，并设置透明度
            highlight_color = [184/255, 160/255, 52/255];%枯萎！！
%             highlight_color = [255/255, 186/255, 58/255];
            p_fill = patch(xverts_fill, yverts_fill, highlight_color, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 无边界颜色

            plot(drought_period, drought_T_rate, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off'); % 底部边界线

            plot(drought_period([1, end]), drought_T_rate([1, end]), 'p', ...
                'MarkerSize', 3, 'LineWidth', 2, 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerFaceColor', [181/255, 18/255, 16/255], 'HandleVisibility', 'off');
            
            % 设置标题和 y 轴标签，根据 third_word 和 drought_speed
            if isempty(tree_kind)
                if strcmp(third_word, 'FD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'NEP')
                        title_str = sprintf('%s Response to %s Flash Drought', data_kind, season_file);
                        ylabel_str = '';
                    end
                elseif strcmp(third_word, 'SD')
                    if strcmp(data_kind, 'T')
                        title_str = sprintf('%s Response to %s Slow Drought', data_kind, season_file);
                        ylabel_str = 'mm day⁻¹';
                    elseif strcmp(data_kind, 'NEP')
                        title_str = sprintf('%s Response to %s Slow Drought', data_kind, season_file);
                        ylabel_str = '';
                    end
                end
            else
                if strcmp(third_word, 'FD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NEP')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NEP')
                            title_str = sprintf('%s Response to %s Flash Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    end
                elseif strcmp(third_word, 'SD')
                    if strcmp(tree_kind, 'Natural')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NEP')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    elseif strcmp(tree_kind, 'Planted')
                        if strcmp(data_kind, 'T')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = 'mm day⁻¹';
                        elseif strcmp(data_kind, 'NEP')
                            title_str = sprintf('%s Response to %s Slow Drought(%s)', data_kind, season_file, tree_kind);
                            ylabel_str = '';
                        end
                    end
                end
            end

            % 获取当前坐标轴对象
            ax = gca;
            set(gca,  'XTickLabel', []);
            
            % 获取坐标轴位置信息
            ax_position = ax.Position;
            
            ax.YAxis.Label.Position(1) = -1.5;
            
            % 设置 x 轴和 y 轴的线条粗细
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            ax.YAxis.FontWeight = 'bold';

                                                % 获取当前子图的当前位置
            pos = get(gca, 'Position');
            
            % 增加高度并保持上端位置不变
            new_height = pos(4) * 1.3;  % 增加高度
            pos(2) = pos(2) - (new_height - pos(4));  % 调整下端位置
            pos(4) = new_height;  % 设置新的高度
            
            % 应用新的位置
            set(gca, 'Position', pos);

            title_handle_1 = title(title_str, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            ylabel(ylabel_str, 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');

            if strcmp(season_file, 'Spring')
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Transpiration','Transpiration Anomaly'}, 'Location', 'northwest','FontSize', 14, 'FontWeight', 'bold');
            else
                legend_handle = legend([p_star, p_drought, p_original, p_no_drought, p_later, p_fill], {'Start and end point', 'Drought', 'Original', 'No Drought', 'Cumulative Transpiration','Transpiration Anomaly'}, 'Location', 'northeast','FontSize', 18, 'FontWeight', 'bold');
            end
            legend_handle.Box = 'off';

            y_norm2 = 0.06/3;
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
            y_norm3 = 0.15/3;
            y_relative3 = y_norm3 * ax_position(4);

            box off     % 取消边框
            ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax1.XAxis.LineWidth = 2;
            ax1.YAxis.LineWidth = 2;
            ax1.YAxis.FontWeight = 'bold';

            

            % 第二个子图：绘制新的时间点区间的蒸腾变化量曲线和填充色
            h2 = subplot(2, 1, 2);
            plot(periods_before_drought, cumulative_transpiration_before_drought, '-k','LineWidth', 2.5);
            hold on
            
            offset = cumulative_transpiration_before_drought(end); % 你想加上的值
            
            % 修改后的数据
            cumulative_transpiration_drought_offset = cumulative_transpiration_droughting + offset;
            cumulative_transpiration_lush_offset = cumulative_transpiration_lush + offset;
            
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5);
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [12/255, 114/255, 31/255],'LineWidth',2.5);

            % 根据 fourth_word 设置不同的 x 和 y 轴范围和刻度间隔
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Spring'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 10];    
                            yticks_custom2 = 0:2:10;
                        case 'Summer'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 50];    
                            yticks_custom2 = 0:10:50;
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 40];    
                            yticks_custom2 = 0:10:40;
                        case 'Winter'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 12];    
                            yticks_custom2 = 0:4:12;
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Spring'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 16];    
                            yticks_custom2 = 0:4:16;
                        case 'Summer'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 50];    
                            yticks_custom2 = 0:10:50;
                        case 'Autumn'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 40];    
                            yticks_custom2 = 0:10:40;
                        case 'Winter'
                            xlim_custom2 = [0 26];        
                            xticks_custom2 = 0:1:26;      
                            ylim_custom2 = [0 10];    
                            yticks_custom2 = 0:2:10;
                    end
                end
            else
                switch season_file
                    case 'Spring'
                        xlim_custom2 = [0 26];
                        xticks_custom2 = 0:1:26;
                        ylim_custom2 = [0 10];
                        yticks_custom2 = 0:2:10;
                    case 'Summer'
                        xlim_custom2 = [0 26];        
                        xticks_custom2 = 0:1:26;      
                        ylim_custom2 = [0 45];    
                        yticks_custom2 = 0:15:45;
                    case 'Autumn'
                        xlim_custom2 = [0 26];        
                        xticks_custom2 = 0:1:26;      
                        ylim_custom2 = [0 40];    
                        yticks_custom2 = 0:10:40;
                    case 'Winter'
                        xlim_custom2 = [0 26];        
                        xticks_custom2 = 0:1:26;      
                        ylim_custom2 = [0 12];    
                        yticks_custom2 = 0:4:12;
                end
            end

            % 设置 x 和 y 轴范围
            xlim(xlim_custom2);
            ylim(ylim_custom2);
            
            % 设置 x 和 y 轴刻度间隔
            xticks(xticks_custom2);
            yticks(yticks_custom2);
            
            % 创建阴影区域
            x_fill_shadow_1 = [period_points(1), serious_period, serious_period, period_points(1)];
            y_fill_shadow_1 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_1, y_fill_shadow_1, [128/255, 128/255, 128/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_2 = [serious_period, period_points(2), period_points(2), serious_period];
            y_fill_shadow_2 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_2, y_fill_shadow_2, [128/255, 128/255, 128/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 创建阴影区域
            x_fill_shadow_3 = [period_points(2), 26, 26, period_points(2)];
            y_fill_shadow_3 = [min(ylim), min(ylim), max(ylim), max(ylim)]; % 使用 y 轴的限制
            
            % 添加阴影
            fill(x_fill_shadow_3, y_fill_shadow_3, [128/255, 128/255, 128/255], 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility', 'off');
            
            % 获取填充区域的顶点坐标
            x_fill = [drought_period, fliplr(lush_period)];
            y_fill = [cumulative_transpiration_drought_offset, fliplr(cumulative_transpiration_lush_offset)];
            
            % 填充区域
            fill(x_fill, y_fill, [184/255, 160/255, 52/255], 'FaceAlpha', 0.8, 'EdgeColor', 'none');

            hold on
            plot(drought_period, cumulative_transpiration_drought_offset, '-', 'Color', [181/255, 18/255, 16/255], 'LineWidth', 2.5, 'HandleVisibility', 'off');
            plot(lush_period, cumulative_transpiration_lush_offset, '-','Color', [12/255, 114/255, 31/255],'LineWidth',2.5, 'HandleVisibility', 'off');
            
            % 添加箭头、短线
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [12/255, 114/255, 31/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [12/255, 114/255, 31/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.6, '^r', ...
                            'MarkerFaceColor', [12/255, 114/255, 31/255], 'MarkerEdgeColor', [12/255, 114/255, 31/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+1, 'vr', ...
                            'MarkerFaceColor', [12/255, 114/255, 31/255], 'MarkerEdgeColor', [12/255, 114/255, 31/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [12/255, 114/255, 31/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-1, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.6, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);
                    end

                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        % 对于不干旱的情况
                        % 绘制横线
                        line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [12/255, 114/255, 31/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [12/255, 114/255, 31/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.7, '^r', ...
                            'MarkerFaceColor', [12/255, 114/255, 31/255], 'MarkerEdgeColor', [12/255, 114/255, 31/255], 'MarkerSize', 5);  % 上箭头
                        plot(lush_period(1), cumulative_transpiration_lush_offset(1)+1, 'vr', ...
                            'MarkerFaceColor', [12/255, 114/255, 31/255], 'MarkerEdgeColor', [12/255, 114/255, 31/255],'MarkerSize', 5);  % 下箭头
                        
                        % 计算线段长度
                        length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                        
                        % 添加长度标注
                        text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                            sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [12/255, 114/255, 31/255]);
                        
                        % 对于干旱情况
                        % 绘制横线
                        line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                        % 绘制竖线
                        line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                        % 添加箭头符号
                        plot(drought_period(end), cumulative_transpiration_drought_offset(end)-1.1, '^r', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                        plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.6, 'vr', ...
                            'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                        
                        
                        % 计算线段长度
                        length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                        
                        % 添加长度标注
                        text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                            sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                            'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    % 对于不干旱的情况
                    % 绘制横线
                    line([min(lush_period), max(lush_period)], [cumulative_transpiration_lush_offset(end), cumulative_transpiration_lush_offset(end)], 'Color', [12/255, 114/255, 31/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([lush_period(1), lush_period(1)], [cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)], 'Color', [12/255, 114/255, 31/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(lush_period(1), cumulative_transpiration_lush_offset(end)-0.7, '^r', ...
                        'MarkerFaceColor', [12/255, 114/255, 31/255], 'MarkerEdgeColor', [12/255, 114/255, 31/255], 'MarkerSize', 5);  % 上箭头
                    plot(lush_period(1), cumulative_transpiration_lush_offset(1)+0.9, 'vr', ...
                        'MarkerFaceColor', [12/255, 114/255, 31/255], 'MarkerEdgeColor', [12/255, 114/255, 31/255],'MarkerSize', 5);  % 下箭头
                    
                    % 计算线段长度
                    length_of_line = abs(cumulative_transpiration_lush_offset(end) - cumulative_transpiration_lush_offset(1));
                    
                    % 添加长度标注
                    text(lush_period(1)-0.25, mean([cumulative_transpiration_lush_offset(1), cumulative_transpiration_lush_offset(end)]), ...
                        sprintf('%.2f', length_of_line), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [12/255, 114/255, 31/255]);
                    
                    % 对于干旱情况
                    % 绘制横线
                    line([min(drought_period), max(drought_period)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(1)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 0.8);
                    % 绘制竖线
                    line([drought_period(end), drought_period(end)], [cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)], 'Color', [181/255, 18/255, 16/255], 'LineWidth', 1.6);
                    % 添加箭头符号
                    plot(drought_period(end), cumulative_transpiration_drought_offset(end)-0.8, '^r', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 上箭头
                    plot(drought_period(end), cumulative_transpiration_drought_offset(1)+0.6, 'vr', ...
                        'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], 'MarkerSize', 5);  % 下箭头
                    
                    
                    % 计算线段长度
                    length_of_line2 = abs(cumulative_transpiration_drought_offset(end) - cumulative_transpiration_drought_offset(1));
                    
                    % 添加长度标注
                    text(drought_period(end)+0.25, mean([cumulative_transpiration_drought_offset(1), cumulative_transpiration_drought_offset(end)]), ...
                        sprintf('%.2f', length_of_line2), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                        'FontName', 'Times New Roman', 'FontSize', 18,'FontWeight', 'bold','Color', [181/255, 18/255, 16/255]);
                end
            end
            
            scatter(drought_period(1), cumulative_transpiration_drought_offset(1), ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(drought_period(end), cumulative_transpiration_drought_offset(end), ...
                'MarkerFaceColor', [181/255, 18/255, 16/255], 'MarkerEdgeColor', [181/255, 18/255, 16/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            scatter(lush_period(end), cumulative_transpiration_lush_offset(end), ...
                'MarkerFaceColor', [12/255, 114/255, 31/255], 'MarkerEdgeColor', [12/255, 114/255, 31/255], ...
                'LineWidth', 1, 'SizeData', 30, 'HandleVisibility', 'off');
            
            % 黑色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.3, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{——}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                            text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.3, ...
                                ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{——}'], ...
                                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                                'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                        text(drought_period(1), cumulative_transpiration_drought_offset(1)+0.3, ...
                            ['{\bf ', num2str(round(cumulative_transpiration_drought_offset(1), 2)), '} \bf{——}'], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                            'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
                end
            end

            % 红色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.04, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                        text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.2, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                    text(drought_period(end), cumulative_transpiration_drought_offset(end)+0.02, ...
                        ['\bf{——} ', num2str(round(cumulative_transpiration_drought_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [181/255, 18/255, 16/255]);
                end
            end

            % 绿色标注
            if ~isempty(tree_kind)
                if strcmp(tree_kind,'Natural')
                    switch season_file
                        case 'Autumn'
                        text(lush_period(end), cumulative_transpiration_lush_offset(end)+1.8, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [12/255, 114/255, 31/255]);
                    end
                elseif strcmp(tree_kind,'Planted')
                    switch season_file
                        case 'Autumn'
                         text(lush_period(end), cumulative_transpiration_lush_offset(end)+1.8, ...
                            ['\bf{—} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                            'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [12/255, 114/255, 31/255]);
                    end
                end
            else
                switch season_file
                    case 'Autumn'
                     text(lush_period(end), cumulative_transpiration_lush_offset(end)+2.1, ...
                        ['\bf{——} ', num2str(round(cumulative_transpiration_lush_offset(end), 2))], ...
                        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontName', ...
                        'Times New Roman','FontSize', 20, 'FontWeight', 'bold', 'Color', [12/255, 114/255, 31/255]);
                end
            end

            
            % 获取当前坐标轴对象
            ax2 = gca;
            set(gca,  'XTickLabel', []);

            % 设置 x 轴和 y 轴的线条粗细
            ax2.XAxis.LineWidth = 2;
            ax2.YAxis.LineWidth = 2;
            ax2.XAxis.FontSize = 20;
            ax2.YAxis.FontSize = 20;
            ax2.YAxis.FontWeight = 'bold';

            title_handle_2 = title('Cumulative Transpiration', 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');
            xlabel('Periods', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
            ylabel('mm', 'FontName', 'Times New Roman','FontSize', 24, 'FontWeight', 'bold');
            
            if strcmp(tree_kind, 'Planted') || strcmp(tree_kind, 'Natural')
                h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 14, 'FontWeight', 'bold');
            else
                if strcmp(season_file,'Autumn')
                    h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northeast', 'FontSize', 18, 'FontWeight', 'bold');
                else
                    h_legend = legend('Original','Drought', 'No Drought', 'Location', 'northwest', 'FontSize', 14, 'FontWeight', 'bold');
                end
            end
            h_legend.Box = 'off';
            
            % 获取坐标轴位置信息
            ax2_position2 = ax2.Position;

                                    % 调整 x 轴标签位置
            y_norm2 = 0.075/3 ;
            y_relative = y_norm2 * ax_position(4);

            abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            % 在原有刻度上添加额外的标记
            for i = 1:length(positions)
                x_position = positions(i);
                x_norm2 = (x_position - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));

                x_relative = x_norm2 * ax2_position2(3);

                % 绘制垂直于 x 轴的线条
                line([x_relative*(26/ax2_position2(3)), x_relative*(26/ax2_position2(3))], [abs_length + ylim_custom2(1),ylim_custom2(1)], ...
                    'Color', 'k', 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
            
            % 调整 x 轴标签位置
            abs_length3 = y_relative3*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));

            ax2.XAxis.Label.Position(2) = ax2.XAxis.Label.Position(2) - abs_length3;

            % 循环添加标签
            for i = 1:length(positions)
                x_label = positions(i);
                x_norm_label = (x_label - xlim_custom2(1)) / (xlim_custom2(2) - xlim_custom2(1));
        
                x_relative_label = x_norm_label * ax2_position2(3);
        
                abs_length = y_relative*((ylim_custom2(2)-ylim_custom2(1))/ax2_position2(4));
                
                % 添加标签
                text(x_relative_label*(26/ax2_position2(3)), ylim_custom2(1)-abs_length, labels{i}, 'Color', 'k', 'VerticalAlignment', 'top', ...
                    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
            end
            
            ax2.YAxis.Label.Position(1) = -1.5;
            
            % 调整垂直间距
            vertical_gap = 0.06;
            % 获取第一个子图的位置
            pos1 = get(h1, 'Position');
            
            % 获取第二个子图的位置，并调整其位置
            pos2 = get(h2, 'Position');
            pos2(2) = pos1(2) - pos1(4) + vertical_gap; % 调整第二个子图的 y 坐标位置
            
            % 更新第二个子图的位置
            set(h2, 'Position', pos2);

            box off     % 取消边框
            ax2_2 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
                'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
            set(ax2_2,'XTick', [],'YTick', []);   % 去掉xy轴刻度

            % 设置 x 轴和 y 轴的线条粗细
            ax2_2.XAxis.LineWidth = 2;
            ax2_2.YAxis.LineWidth = 2;
            ax2_2.YAxis.FontWeight = 'bold';
            
            % 设置保存路径和文件名
            save_path = ['E:\Data\Drought_quantification\Quantify(pictures)\Quantify_0403\' segmented_name];
            
            % 完整的保存路径
            full_file_path = fullfile(save_path, output_file);
            
            % 设置输出分辨率为 300 DPI
            resolution = 500;
            
            % 保存图形为 JPEG 格式
            saveas(gcf, full_file_path, 'svg');
        end
    end
end
