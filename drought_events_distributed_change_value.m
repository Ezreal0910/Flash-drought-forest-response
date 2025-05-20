clc;
clear;

% 读取栅格数据
[A,R]=readgeoraster('F:\GZW\Output\Replace.tif');

period = 1:25;
serious_period = 13;

% 定义新的输出子文件夹名称
output_path = 'L:\GZW\Drought_quantify\Drought_distributed2_difference';

% 定义输入文件夹路径列表
input_folders = {
    'L:\GZW\Drought_quantify\Season_ALL_1_month\FD',...
    'L:\GZW\Drought_quantify\Season_N&P_1_month\Natural\FD',...
    'L:\GZW\Drought_quantify\Season_N&P_1_month\Planted\FD'
};

% 循环处理每个输入文件夹
for f = 1:length(input_folders)
    folder_path = input_folders{f};
    input_path = fullfile(folder_path, '1CI_new');
    
    % 定义四个季节的名称
    seasons = {'spring', 'summer', 'autumn', 'winter'};
    
    % 初始化一个 cell 数组来存储季节文件夹的完整路径
    season_paths = cell(size(seasons));
    
    % 构建每个季节文件夹的完整路径
    for i = 1:numel(seasons)
        season_paths{i} = fullfile(input_path, seasons{i});
    end
    
    % 现在您可以在每个季节文件夹路径上运行之前修改过的处理代码
    for season_idx = 1:numel(season_paths) 
        season_folder = season_paths{season_idx};
    
        [~, season_file] = fileparts(season_folder);
    
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
                % 获取子文件夹名称
                [~, folder_name_short, ~] = fileparts(folder_name);
    
                % 检查子文件夹是否为 'Mean'
                if strcmp(folder_name_short, 'Mean')
                    disp(['读取处理 文 件 夹 ：', folder_name]);
                
                    % 获取当前文件夹中的所有文件和文件夹
                    sub_contents = dir(folder_name);
    
                    % 初始化 cell 数组来存储文件名
                    selected_files_T = {};
                    selected_files_T_start_periods = {};
                    selected_files_T_end_periods = {};
    
                    selected_files_NEP = {};
                    selected_files_NEP_start_periods = {};
                    selected_files_NEP_end_periods = {};
    
                    selected_files_GPP = {};
                    selected_files_GPP_start_periods = {};
                    selected_files_GPP_end_periods = {};

                    selected_files_NPP = {};
                    selected_files_NPP_start_periods = {};
                    selected_files_NPP_end_periods = {};

                    selected_files_Ra = {};
                    selected_files_Ra_start_periods = {};
                    selected_files_Ra_end_periods = {};
    
                    selected_files_Reco = {};
                    selected_files_Reco_start_periods = {};
                    selected_files_Reco_end_periods = {};
            
                    % 遍历当前文件夹中的所有内容
                    for j = 1:numel(sub_contents)
                        % 忽略当前目录（'.'）和父目录（'..'）
                        if strcmp(sub_contents(j).name, '.') || strcmp(sub_contents(j).name, '..')
                            continue;
                        end
    
                        file_name = sub_contents(j).name;
                        file_parts = split(file_name, '_');  % 假设文件名以 '_' 分隔各个部分
    
                        fifth_word = file_parts{4};
                        sixth_word = file_parts{5};
    
                        % 将文件名归类到 cell 数组中
    
                        if strcmp(fifth_word, 'T')
                            if contains(sixth_word, 'start')
                                selected_files_T_start_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            elseif contains(sixth_word, 'end')
                                selected_files_T_end_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            else
                                selected_files_T{end+1} = fullfile(folder_name, sub_contents(j).name);
                            end
                        elseif strcmp(fifth_word, 'NEP')
                            if contains(sixth_word, 'start')
                                selected_files_NEP_start_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            elseif contains(sixth_word, 'end')
                                selected_files_NEP_end_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            else
                                selected_files_NEP{end+1} = fullfile(folder_name, sub_contents(j).name);
                            end
                        elseif strcmp(fifth_word, 'GPP')
                            if contains(sixth_word, 'start')
                                selected_files_GPP_start_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            elseif contains(sixth_word, 'end')
                                selected_files_GPP_end_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            else
                                selected_files_GPP{end+1} = fullfile(folder_name, sub_contents(j).name);
                            end
                        elseif strcmp(fifth_word, 'NPP')
                            if contains(sixth_word, 'start')
                                selected_files_NPP_start_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            elseif contains(sixth_word, 'end')
                                selected_files_NPP_end_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            else
                                selected_files_NPP{end+1} = fullfile(folder_name, sub_contents(j).name);
                            end
                        elseif strcmp(fifth_word, 'Ra')
                            if contains(sixth_word, 'start')
                                selected_files_Ra_start_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            elseif contains(sixth_word, 'end')
                                selected_files_Ra_end_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            else
                                selected_files_Ra{end+1} = fullfile(folder_name, sub_contents(j).name);
                            end
                        elseif strcmp(fifth_word, 'Reco')
                            if contains(sixth_word, 'start')
                                selected_files_Reco_start_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            elseif contains(sixth_word, 'end')
                                selected_files_Reco_end_periods{end+1} = fullfile(folder_name, sub_contents(j).name);
                            else
                                selected_files_Reco{end+1} = fullfile(folder_name, sub_contents(j).name);
                            end
                        end
                    end
                end
            end
        end

        % 将这些变量放入一个 cell 数组中
        selected_files_all = {selected_files_T, selected_files_NEP, selected_files_GPP, selected_files_NPP, selected_files_Ra, selected_files_Reco};
        selected_files_all_start = {selected_files_T_start_periods, selected_files_NEP_start_periods, selected_files_GPP_start_periods, selected_files_NPP_start_periods, selected_files_Ra_start_periods, selected_files_Reco_start_periods};
        selected_files_all_end = {selected_files_T_end_periods, selected_files_NEP_end_periods, selected_files_GPP_end_periods, selected_files_NPP_end_periods, selected_files_Ra_end_periods, selected_files_Reco_end_periods};
    
        for m = 1:numel(selected_files_all)
            selected_file = selected_files_all{m};
            selected_file = selected_file{1};

            selected_file_start = selected_files_all_start{m};
            selected_file_start = selected_file_start{1};

            selected_file_end = selected_files_all_end{m};
            selected_file_end = selected_file_end{1};
            
            [~, tree_kind] = fileparts(fileparts(fileparts(fileparts(fileparts(fileparts(selected_file))))));
    
            if contains(tree_kind, 'Natural')
                tree_kind = 'Natural';
            elseif contains(tree_kind, 'Planted')
                tree_kind = 'Planted';
            elseif contains(tree_kind, 'ALL')
                tree_kind = 'ALL';
            end
    
            % 提取文件名部分（不包括路径和扩展名）
            [~, data_kinds, ~] = fileparts(selected_file);
            data_words = split(data_kinds, '_');
            % 获取第四个单词
            third_word = data_words{3};
            data_kind = data_words{4};

            selected_file_path_parts = split(selected_file, filesep);
            fourth_file_name = selected_file_path_parts(4);

            savePath1 = fullfile(output_path, tree_kind, third_word, season_file);

            % 检查保存路径是否存在，不存在则创建
            if ~isfolder(savePath1)
                mkdir(savePath1); % 创建保存路径
            end

            if strcmp(fourth_file_name, 'Season_ALL_1_month')
                % 构造对应文件的路径
                corresponding_file = strrep(selected_file, 'Season_ALL_1_month', 'Season_ALL_no_drought_1_month');
                path_parts = split(corresponding_file, filesep);
                path_parts(8) = [];

                output_file = sprintf('%s_%s_%s_%s_difference.tif', tree_kind, data_kind, third_word, season_file);
            else
                corresponding_file = strrep(selected_file, 'Season_N&P_1_month', 'Season_N&P_no_drought_1_month');
                path_parts = split(corresponding_file, filesep);
                path_parts(9) = [];

                % 构建新的输出文件名
                output_file = sprintf('%s_%s_%s_%s_difference.tif', tree_kind, data_kind, third_word, season_file);
            end
            % 重新组合路径
            corresponding_file = fullfile(path_parts{:});
    
            data_T = load(selected_file);
            data_T_names = fieldnames(data_T);
            data_T_name = data_T_names{1};
            T_data_3D = data_T.(data_T_name); 

            data_T_no_drought = load(corresponding_file);
            data_T_names_no_drought = fieldnames(data_T_no_drought);
            data_T_name_no_drought = data_T_names_no_drought{1};
            T_data_3D_no_drought = data_T_no_drought.(data_T_name_no_drought);

            data_T_start = load(selected_file_start);
            data_T_start_names = fieldnames(data_T_start);
            data_T_start_name = data_T_start_names{1};
            T_data_start_2D = data_T_start.(data_T_start_name);

            data_T_end = load(selected_file_end);
            data_T_end_names = fieldnames(data_T_end);
            data_T_end_name = data_T_end_names{1};
            T_data_end_2D = data_T_end.(data_T_end_name);

            forest_FD_T_Difference_value = zeros(3900, 6200);

            % 找到非零值的位置
            [rows2, cols2] = find(T_data_start_2D);

            for i = 1:length(rows2)
                row2 = rows2(i);
                col2 = cols2(i);


                data_T_values = T_data_3D(row2, col2, :);
                data_T_values = data_T_values(:);
                data_T_values = data_T_values(:)';

                if any(isnan(data_T_values))  
                    continue;  
                end 
            
                drought_start = T_data_start_2D(row2, col2);
                drought_end = T_data_end_2D(row2, col2);

                % 使用 interp1 函数进行插值
                value_start = interp1(period, data_T_values, drought_start, 'linear');
                value_end = interp1(period, data_T_values, drought_end, 'linear');

                if isnan(value_start) || isnan(value_end)  
                    continue;
                end

                data_T_values_no_drought = T_data_3D_no_drought(row2, col2, :);
                data_T_values_no_drought = data_T_values_no_drought(:);
                data_T_values_no_drought = data_T_values_no_drought(:)';

                if any(isnan(data_T_values_no_drought))  
                    continue;  
                end 

                % 使用 interp1 函数进行插值
                value_start_no_drought = interp1(period, data_T_values_no_drought, drought_start, 'linear');
                value_end_no_drought = interp1(period, data_T_values_no_drought, drought_end, 'linear');

                if isnan(value_start_no_drought) || isnan(value_end_no_drought)  
                    continue;
                end
    
                % 提取从 start_period 开始的时间和蒸腾变化量数据
                idx_start = find(period >= ceil(drought_start), 1);
                idx_end = find(period >= ceil(drought_end), 1);
                drought_period = [drought_start, period(idx_start:idx_end-1),drought_end];
                drought_T_rate = [value_start, data_T_values(idx_start:idx_end-1),value_end];
                cumulative_transpiration_droughting = cumtrapz(drought_period, drought_T_rate);

                integer_periods = ceil(drought_start):floor(drought_end);
                
                lush_period = [drought_start, integer_periods, drought_end];
                lush_T_rate = [value_start_no_drought, data_T_values_no_drought(integer_periods), value_end_no_drought]; 
                cumulative_transpiration_lush = cumtrapz(lush_period, lush_T_rate);

                Difference_value = cumulative_transpiration_lush(end) - cumulative_transpiration_droughting(end);

                forest_FD_T_Difference_value(row2, col2) = Difference_value;
            end

            % 目标后缀
            csv_extension = '.csv';
            mat_extension = '.mat';
            
            % 找到最后一个点的位置
            dot_index = strfind(output_file, '.');
            
            % 如果存在点，则替换后缀
            if ~isempty(dot_index)
                % 创建新的文件名
                csv_file_name = [output_file(1:dot_index(end)-1) csv_extension];
                mat_file_name = [output_file(1:dot_index(end)-1) mat_extension];
            else
                % 如果没有点，则直接添加新的后缀
                csv_file_name = [output_file csv_extension];
                mat_file_name = [output_file mat_extension];
            end

            save(fullfile(savePath1, mat_file_name), 'forest_FD_T_Difference_value', '-v7.3');

            % 保存结果为CSV文件
            num_more_than_zero = sum(forest_FD_T_Difference_value(:) > 0);
            num_non_zero = sum(forest_FD_T_Difference_value(:) ~= 0);
            
            % 计算大于0的数值占不为零数值的比例
            proportion_more_than_zero = num_more_than_zero / num_non_zero;
            count_results = [num_more_than_zero, num_non_zero, proportion_more_than_zero];
            writematrix(count_results, fullfile(savePath1, csv_file_name));

            % 写回栅格数据
            new_output_path = fullfile(savePath1, output_file);
            geotiffwrite(new_output_path, forest_FD_T_Difference_value, R);
            fprintf('已输出 %s : %s_forest_%s_%s\n', season_file, tree_kind, third_word, data_kind);
            
            % 输出结果
            fprintf('大于0的数值个数: %d\n', num_more_than_zero);
            fprintf('大于0的数值占不为零数值的比例: %.2f\n', proportion_more_than_zero);
        end
    end
end