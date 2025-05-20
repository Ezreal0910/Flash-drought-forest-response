clc;
clear;

% 要计算的置信水平数组
confidence_levels = 5:10:95;

% 转换为累积概率
cumulative_probabilities = (100 + confidence_levels) / 200;

% 计算对应的 z 值
z_values = norminv(cumulative_probabilities);

% 定义输入文件夹路径列表
input_folders = {
    'M:\GZW\Drought_quantify\Season_ALL_1_month\FD\1CI_new',...
    'M:\GZW\Drought_quantify\Season_N&P_1_month\Natural\FD\1CI_new',...
    'M:\GZW\Drought_quantify\Season_N&P_1_month\Planted\FD\1CI_new'
};

% 循环处理每个输入文件夹
for f = 1:length(input_folders)
    folder_path = input_folders{f};
    output_path = fullfile(folder_path, '1CI');

    % 获取文件夹下所有子文件夹名称
    subfolders = dir(folder_path);
    subfolders = subfolders([subfolders.isdir]);  % 仅保留文件夹
    
    % 用于存储子文件夹名称的单元数组
    folder_names = {};
    
    % 遍历子文件夹并提取名称
    for i = 1:length(subfolders)
        % 排除当前目录（.）和上级目录（..）
        if ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..')
            folder_names{end+1} = subfolders(i).name;
        end
    end
    
    for j = 1:length(folder_names)
        current_folder_name = folder_names{j};
        if strcmp(current_folder_name, '1CI')
            continue;
        end
        current_folder_path = fullfile(folder_path, current_folder_name);

        files = dir(fullfile(current_folder_path, '*'));
    
        % 存储文件完整路径的数组
        file_paths = {};
        
        for m = 1:numel(files)
            % 检查文件名是否不是 "." 和 ".."，且不含有数字，并且是文件而不是文件夹
            if ~strcmp(files(m).name, '.') && ~strcmp(files(m).name, '..') && ...
                    isempty(regexp(files(m).name, '\d', 'once')) && ...
                    files(m).isdir
                % 构建完整路径并存储到数组
                file_paths{end+1} = fullfile(current_folder_path, files(m).name);
            end
        end

        % 遍历符合条件的文件夹
        for n = 1:length(file_paths)
            current_folder_path = file_paths{n};
            disp(['Processing folder: ' current_folder_path]);

            [ ~, drought_speed] = fileparts(current_folder_path);
            
            % 获取当前文件夹中的所有文件
            subfiles = dir(fullfile(current_folder_path, '*.mat'));  % 仅获取 .mat 文件

            % 创建一个结构体数组来存储分类后的文件
            data_kinds = {'T', 'NEP', 'GPP', 'NPP', 'Ra', 'Reco'};
            classified_files = struct();
            
            for q = 1:length(data_kinds)
                classified_files.(data_kinds{q}) = {};
            end
            
            
            for r = 1:length(subfiles)
                filename = subfiles(r).name;
                filename_words = split(filename,'_');
        
                if length(filename_words) >= 4
                    kind_name = filename_words(4);
                    kind_name = kind_name{1};
            
                    % 检查提取出的单词是否为月份，并分类
                    if ismember(kind_name, data_kinds)
                        classified_files.(kind_name){end+1} = filename;
                    else
                        disp(['File: ' filename ', Second word is not a month']);
                    end
                else
                    disp(['File: ' filename ', No second word found']);
                end
            end
            
            % 遍历分类后的文件
            for s = 1:length(data_kinds)
                data_kind = data_kinds{s};
                data_kind_files = classified_files.(data_kind);
       
                filepath_data = fullfile(current_folder_path,data_kind_files{1});
                
                % 加载文件
                loaded_data = load(filepath_data);
                data_names = fieldnames(loaded_data);
                data_variable = loaded_data.(data_names{1});

                % 如果是 double 类型的数据，则执行另一种处理
                disp(['正在处理文件： ', data_kind_files{1}]);
                if strcmp(drought_speed, 'Slow')
                    savePath1 = fullfile(output_path, current_folder_name, 'Slow');
                elseif strcmp(drought_speed, 'Fast')
                    savePath1 = fullfile(output_path, current_folder_name, 'Fast');
                elseif strcmp(drought_speed, 'Mean')
                    savePath1 = fullfile(output_path, current_folder_name, 'Mean');
                end


                % 检查保存路径是否存在，不存在则创建
                if ~isfolder(savePath1)
                    mkdir(savePath1); % 创建保存路径
                end
                
                filepath_end = fullfile(current_folder_path,data_kind_files{2});
                periods_data1 = load(filepath_end);
                periods_data_names1 = fieldnames(periods_data1);
                forest_FD_drought_end = periods_data1.(periods_data_names1{1});
        
                filepath_start = fullfile(current_folder_path,data_kind_files{3});
                periods_data2 = load(filepath_start);
                periods_data_names2 = fieldnames(periods_data2);
                forest_FD_drought_start = periods_data2.(periods_data_names2{1});

                data_variable(data_variable == 0) = nan;
                forest_FD_drought_start(forest_FD_drought_start == 0) = nan;
                forest_FD_drought_end(forest_FD_drought_end == 0) = nan;
                
                forest_FD_drought_start_mean = nanmean(forest_FD_drought_start, 'all');
                forest_FD_drought_end_mean = nanmean(forest_FD_drought_end, 'all');
                
                drought_periods_means = [forest_FD_drought_start_mean; forest_FD_drought_end_mean];
                
                % 预分配数组用于存储结果
                forest_FD_nanmean_1D = zeros(25, 1);
                forest_FD_nanstd_1D = zeros(25, 1);
                
                forest_FD_standard_error_1D = zeros(25, 1);
                forest_FD_count = zeros(25, 1);
                
                % 预分配数组用于存储结果
                forest_FD_confidence_intervals_upper = zeros(25, length(confidence_levels));
                forest_FD_confidence_intervals_lower = zeros(25, length(confidence_levels));
                
                for slice_index = 1:25
                    % 逐个提取切片
                    current_slice = data_variable(:, :, slice_index);
                
                    % 在切片上计算均值
                    forest_FD_nanmean_1D(slice_index) = nanmean(current_slice, 'all');
                
                    % 在切片上计算标准差
                    forest_FD_nanstd_1D(slice_index) = nanstd(current_slice, 0, 'all');
                
                    % 在切片上计算标准误差
                    num_observation = sum(~isnan(current_slice(:)));  % 计算非 NaN 值的观测数量
                    forest_FD_count(slice_index) = num_observation;
                
                    forest_FD_standard_error_1D(slice_index) = forest_FD_nanstd_1D(slice_index) / sqrt(num_observation);
                
                    % 计算不同置信度下的置信区间
                    for level_index = 1:length(confidence_levels)
                        % 计算置信区间的上限和下限
                        upper_limit = forest_FD_nanmean_1D(slice_index) + z_values(level_index) * forest_FD_standard_error_1D(slice_index);
                        lower_limit = forest_FD_nanmean_1D(slice_index) - z_values(level_index) * forest_FD_standard_error_1D(slice_index);
                        
                        forest_FD_confidence_intervals_upper(slice_index, level_index) = upper_limit;
                        forest_FD_confidence_intervals_lower(slice_index, level_index) = lower_limit;
                    end
                end
    
                filename1 = [current_folder_name '_forest_FD_' data_kind '_mean.csv'];
                writematrix(forest_FD_nanmean_1D, fullfile(savePath1, filename1));
                
                filename2 = [current_folder_name '_forest_FD_' data_kind '_ci_upper.csv'];
                writematrix(forest_FD_confidence_intervals_upper, fullfile(savePath1, filename2));
                
                filename3 = [current_folder_name '_forest_FD_' data_kind '_ci_lower.csv'];
                writematrix(forest_FD_confidence_intervals_lower, fullfile(savePath1, filename3));
                
                filename4 = [current_folder_name '_forest_FD_' data_kind '_standard_error.csv'];
                writematrix(forest_FD_standard_error_1D, fullfile(savePath1, filename4));
                
                filename5 = [current_folder_name '_forest_FD_' data_kind '_std.csv'];
                writematrix(forest_FD_nanstd_1D, fullfile(savePath1, filename5));
                
                filename6 = [current_folder_name '_forest_FD_' data_kind '_count.csv'];
                writematrix(forest_FD_count, fullfile(savePath1, filename6));
    
                filename7 = [current_folder_name '_forest_FD_' data_kind '_periods_means.csv'];
                writematrix(drought_periods_means, fullfile(savePath1, filename7));
                
                fprintf('已输出 forest_FD_%s \n', data_kind);
          
            end
        end
    end
end