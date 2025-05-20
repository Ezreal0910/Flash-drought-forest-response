% clc;
% clear;
% 
% % 要计算的置信水平数组
% confidence_levels = 5:10:95;
% 
% % 转换为累积概率
% cumulative_probabilities = (100 + confidence_levels) / 200;
% 
% % 计算对应的 z 值
% z_values = norminv(cumulative_probabilities);
% 
% load("G:\GZW\Forest\Forest_ESI.mat")

% 定义输入文件夹路径列表
input_folders = {
    'M:\GZW\Drought_quantify\Season_ALL_1_month\FD',...
    'M:\GZW\Drought_quantify\Season_N&P_1_month\Natural\FD',...
    'M:\GZW\Drought_quantify\Season_N&P_1_month\Planted\FD'
};

% 循环处理每个输入文件夹
for f = 1:length(input_folders)
    folder_path = input_folders{f};
    output_path = fullfile(folder_path, '1CI_new');

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
        current_folder_path = fullfile(folder_path, current_folder_name);
        
        if strcmp(current_folder_name, 'spring')
            % 在这里执行 "spring" 文件夹的处理操作
            files = dir(fullfile(current_folder_path, '*'));
        
            % 存储文件完整路径的数组
            file_paths = {};
            
            for m = 1:numel(files)
                % 检查文件名是否不是 "." 和 ".."，且不含有数字，并且是文件而不是文件夹
                if ~strcmp(files(m).name, '.') && ~strcmp(files(m).name, '..') && ...
                        isempty(regexp(files(m).name, '\d', 'once')) && ...
                        ~files(m).isdir
                    % 构建完整路径并存储到数组
                    file_paths{end+1} = fullfile(current_folder_path, files(m).name);
                end
            end
    
            for x = 1:numel(file_paths)
                % 加载文件
                loaded_data = load(file_paths{x});
                data_names = fieldnames(loaded_data);
                data_variable = loaded_data.(data_names{1});
                
                % 判断加载的数据类型
                if iscell(data_variable)
                   % 如果是cell数组类型的数据，则执行一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 cell 类型数据 ']);
    
                    periods_folder_path1 = fullfile(folder_path, current_folder_name);
        
                    % 获取 "spring" 文件夹下的文件列表
                    spring_periods_folder_path1 = dir(fullfile(periods_folder_path1, '*'));
                    spring_periods_folder_path2 = fullfile(periods_folder_path1, spring_periods_folder_path1(3).name);
    
                    spring_periods_files = dir(fullfile(spring_periods_folder_path2, '*'));
    
                    
                    % 存储文件完整路径的数组
                    spring_periods_file_paths = {};
        
                    for m = 1:numel(spring_periods_files)
                        % 检查文件名是否不是 "." 和 ".."，并且不包含数字
                        if ~strcmp(spring_periods_files(m).name, '.') && ...
                           ~strcmp(spring_periods_files(m).name, '..') && ...
                           ~any(isstrprop(spring_periods_files(m).name, 'digit'))
                       
                            % 构建完整路径并添加到数组中
                            spring_file_path = fullfile(spring_periods_folder_path2, spring_periods_files(m).name);
                            spring_periods_file_paths{end+1} = spring_file_path;
                        end
                    end
    
                    matched_files = {}; % 用于存储匹配的文件路径
                    % 使用 data_names{1} 构建正则表达式模式
                    pattern = ['^', regexprep(data_names{1}, '_', '\_'), '_\w+$'];
                    
                    for i = 1:numel(spring_periods_file_paths)
                        [~, filename, ~] = fileparts(spring_periods_file_paths{i}); % 获取文件名
                    
                        % 使用正则表达式模式匹配文件名
                        if ~isempty(regexp(filename, pattern, 'once'))
                            matched_files{end+1} = spring_periods_file_paths{i}; % 将匹配的文件路径添加到 matched_files 数组中
                        end
                    end
    
                    savePath1 = fullfile(output_path, 'spring', 'Mean');
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
     
                    periods_data1 = load(matched_files{1});
                    periods_data_names1 = fieldnames(periods_data1);
                    periods_end = periods_data1.(periods_data_names1{1});
    
                    periods_data2 = load(matched_files{2});
                    periods_data_names2 = fieldnames(periods_data2);
                    periods_start = periods_data2.(periods_data_names2{1});
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    forest_FD_ESI_start = zeros(3900, 6200);
                    forest_FD_ESI_end = zeros(3900, 6200);
                    
                    % 创建一个逻辑索引，表示每个单元格是否非空
                    nonEmptyCells = ~cellfun('isempty', data_variable);
                    
                    % 使用 find 函数找到非空元素的位置
                    [rows1, cols1] = find(nonEmptyCells);
                    
                    skippedCount = 0;
                    
                    for i = 1:length(rows1)
                        row1 = rows1(i);
                        col1 = cols1(i);
                        values = data_variable{row1, col1};
                        drought_starts = periods_start{row1, col1};
                        drought_ends = periods_end{row1, col1};
                    
                        % 获取 values 中的值的数量
                        num_values = numel(values);
                        
                        % 创建 25*values数量的矩阵
                        matrix_size = [25, num_values];
                        created_matrix = zeros(matrix_size);
                    
                        FD_drought_start = [];
                        FD_drought_end = [];
                    
                        % 循环处理 values 中的每个值
                        for k = 1:length(values)
                            current_value = values(k);
                    
                            drought_start = drought_starts(k);
                            drought_end = drought_ends(k);
                        
                            drought_develop1 = current_value - drought_start;
                            drought_recover1 = drought_end - current_value;
                        
                            drought_develop2 = 13 - drought_develop1;
                            drought_recover2 = 13 + drought_recover1;
                    
                            FD_drought_start = [FD_drought_start, drought_develop2];
                            FD_drought_end = [FD_drought_end, drought_recover2];
                    
                            % 检查value是否在指定范围内
                            if current_value <= 12 || current_value > 908
                                skippedCount = skippedCount + 1;
                                continue;
                            end
                        
                            % 获得前后8个索引号
                            startIndex =  current_value - 12;
                            endIndex   =  current_value + 12;
                        
                            forest_ESI = forest_ESI_data(row1, col1, startIndex:endIndex);
                    
                            forest_ESI = forest_ESI(:);
                    
                            created_matrix(1:25, k) = forest_ESI;
                        end
                    
                        created_matrix(created_matrix == 0) = nan;
                        % 计算每一行的均值
                        row_means = nanmean(created_matrix, 2);
                        FD_drought_start_mean = nanmean(FD_drought_start);
                        FD_drought_end_mean = nanmean(FD_drought_end);
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row1, col1, 1:25) = row_means;
                        forest_FD_ESI_start(row1, col1) = FD_drought_start_mean;
                        forest_FD_ESI_end(row1, col1) = FD_drought_end_mean;
                    end
                elseif isnumeric(data_variable)
                    % 如果是 double 类型的数据，则执行另一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 double 类型数据 ']);
    
                    periods_folder_path1 = fullfile(folder_path, current_folder_name);
        
                    % 获取 "spring" 文件夹下的文件列表
                    spring_periods_folder_path1 = dir(fullfile(periods_folder_path1, '*'));
                    spring_periods_folder_path2 = fullfile(periods_folder_path1, spring_periods_folder_path1(3).name);
    
                    spring_periods_files = dir(fullfile(spring_periods_folder_path2, '*'));
                    
                    % 存储文件完整路径的数组
                    spring_periods_file_paths = {};
        
                    for m = 1:numel(spring_periods_files)
                        % 检查文件名是否不是 "." 和 ".."，并且不包含数字
                        if ~strcmp(spring_periods_files(m).name, '.') && ...
                           ~strcmp(spring_periods_files(m).name, '..') && ...
                           ~any(isstrprop(spring_periods_files(m).name, 'digit'))
                       
                            % 构建完整路径并添加到数组中
                            spring_file_path = fullfile(spring_periods_folder_path2, spring_periods_files(m).name);
                            spring_periods_file_paths{end+1} = spring_file_path;
                        end
                    end
    
                    matched_files = {}; % 用于存储匹配的文件路径
    
                    % 使用 data_names{1} 构建正则表达式模式
                    pattern = ['^', regexprep(data_names{1}, '_', '\_'), '_\w+$'];
                    
                    for i = 1:numel(spring_periods_file_paths)
                        [~, filename, ~] = fileparts(spring_periods_file_paths{i}); % 获取文件名
                    
                        % 使用正则表达式模式匹配文件名
                        if ~isempty(regexp(filename, pattern, 'once'))
                            matched_files{end+1} = spring_periods_file_paths{i}; % 将匹配的文件路径添加到 matched_files 数组中
                        end
                    end
                    
                    if contains(data_names{1}, 'slow')
                        savePath1 = fullfile(output_path, 'spring', 'Slow');
                    else
                        savePath1 = fullfile(output_path, 'spring', 'Fast');
                    end
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
                    
                    periods_data1 = load(matched_files{1});
                    periods_data_names1 = fieldnames(periods_data1);
                    periods_end = periods_data1.(periods_data_names1{1});
    
                    periods_data2 = load(matched_files{2});
                    periods_data_names2 = fieldnames(periods_data2);
                    periods_start = periods_data2.(periods_data_names2{1});
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    forest_FD_ESI_start = zeros(3900, 6200);
                    forest_FD_ESI_end = zeros(3900, 6200);
    
                    % 找到非零值的位置
                    [rows2, cols2] = find(data_variable);
    
                    for i = 1:length(rows2)
                        row2 = rows2(i);
                        col2 = cols2(i);
                        value = data_variable(row2, col2);
                    
                        drought_start = periods_start(row2, col2);
                        drought_end = periods_end(row2, col2);
                    
                        drought_develop1 = value - drought_start;
                        drought_recover1 = drought_end - value;
                    
                        drought_develop2 = 13 - drought_develop1;
                        drought_recover2 = 13 + drought_recover1;
                    
                        % 检查value是否在指定范围内
                        if value <= 12 || value > 908
                            skippedCount = skippedCount + 1;
                            continue;
                        end
                    
                        % 获得前后8个索引号
                        startIndex =  value - 12;
                        endIndex =  value + 12;
                    
                        forest_ESI = forest_ESI_data(row2, col2, startIndex:endIndex);
                        % 将子集转换为列向量
                        forest_ESI = forest_ESI(:);
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row2, col2, 1:25) = forest_ESI;
                    
                        forest_FD_ESI_start(row2, col2) = drought_develop2;
                        forest_FD_ESI_end(row2, col2) = drought_recover2;
                    end         
                else
                    % 其他数据类型的处理
                    disp(['文件 ', file_paths{x}, ' 不处理']);
                end
    
                filename1 = [current_folder_name '_forest_FD_ESI_data.mat'];
                save(fullfile(savePath1, filename1), 'forest_FD_ESI_data', '-v7.3');

                filename2 = [current_folder_name '_forest_FD_ESI_start.mat'];
                save(fullfile(savePath1, filename2), 'forest_FD_ESI_start', '-v7.3');

                filename3 = [current_folder_name '_forest_FD_ESI_end.mat'];
                save(fullfile(savePath1, filename3), 'forest_FD_ESI_end', '-v7.3');
                
                fprintf('已输出 %s : forest_FD_ESI\n', current_folder_name);
            end
    
        elseif strcmp(current_folder_name, 'summer')
                    % 在这里执行 "summer" 文件夹的处理操作
            files = dir(fullfile(current_folder_path, '*'));
        
            % 存储文件完整路径的数组
            file_paths = {};
            
            for m = 1:numel(files)
                % 检查文件名是否不是 "." 和 ".."，且不含有数字，并且是文件而不是文件夹
                if ~strcmp(files(m).name, '.') && ~strcmp(files(m).name, '..') && ...
                        isempty(regexp(files(m).name, '\d', 'once')) && ...
                        ~files(m).isdir
                    % 构建完整路径并存储到数组
                    file_paths{end+1} = fullfile(current_folder_path, files(m).name);
                end
            end
    
            for x = 1:numel(file_paths)
                % 加载文件
                loaded_data = load(file_paths{x});
                data_names = fieldnames(loaded_data);
                data_variable = loaded_data.(data_names{1});
                
                % 判断加载的数据类型
                if iscell(data_variable)
                   % 如果是cell数组类型的数据，则执行一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 cell 类型数据 ']);
    
                    periods_folder_path1 = fullfile(folder_path, current_folder_name);
        
                    % 获取 "summer" 文件夹下的文件列表
                    summer_periods_folder_path1 = dir(fullfile(periods_folder_path1, '*'));
                    summer_periods_folder_path2 = fullfile(periods_folder_path1, summer_periods_folder_path1(3).name);
    
                    summer_periods_files = dir(fullfile(summer_periods_folder_path2, '*'));
    
                    
                    % 存储文件完整路径的数组
                    summer_periods_file_paths = {};
        
                    for m = 1:numel(summer_periods_files)
                        % 检查文件名是否不是 "." 和 ".."，并且不包含数字
                        if ~strcmp(summer_periods_files(m).name, '.') && ...
                           ~strcmp(summer_periods_files(m).name, '..') && ...
                           ~any(isstrprop(summer_periods_files(m).name, 'digit'))
                       
                            % 构建完整路径并添加到数组中
                            summer_file_path = fullfile(summer_periods_folder_path2, summer_periods_files(m).name);
                            summer_periods_file_paths{end+1} = summer_file_path;
                        end
                    end
    
                    matched_files = {}; % 用于存储匹配的文件路径
                    % 使用 data_names{1} 构建正则表达式模式
                    pattern = ['^', regexprep(data_names{1}, '_', '\_'), '_\w+$'];
                    
                    for i = 1:numel(summer_periods_file_paths)
                        [~, filename, ~] = fileparts(summer_periods_file_paths{i}); % 获取文件名
                    
                        % 使用正则表达式模式匹配文件名
                        if ~isempty(regexp(filename, pattern, 'once'))
                            matched_files{end+1} = summer_periods_file_paths{i}; % 将匹配的文件路径添加到 matched_files 数组中
                        end
                    end
    
                    savePath1 = fullfile(output_path, 'summer', 'Mean');
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
     
                    periods_data1 = load(matched_files{1});
                    periods_data_names1 = fieldnames(periods_data1);
                    periods_end = periods_data1.(periods_data_names1{1});
    
                    periods_data2 = load(matched_files{2});
                    periods_data_names2 = fieldnames(periods_data2);
                    periods_start = periods_data2.(periods_data_names2{1});
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    forest_FD_ESI_start = zeros(3900, 6200);
                    forest_FD_ESI_end = zeros(3900, 6200);
                    
                    % 创建一个逻辑索引，表示每个单元格是否非空
                    nonEmptyCells = ~cellfun('isempty', data_variable);
                    
                    % 使用 find 函数找到非空元素的位置
                    [rows1, cols1] = find(nonEmptyCells);
                    
                    skippedCount = 0;
                    
                    for i = 1:length(rows1)
                        row1 = rows1(i);
                        col1 = cols1(i);
                        values = data_variable{row1, col1};
                        drought_starts = periods_start{row1, col1};
                        drought_ends = periods_end{row1, col1};
                    
                        % 获取 values 中的值的数量
                        num_values = numel(values);
                        
                        % 创建 25*values数量的矩阵
                        matrix_size = [25, num_values];
                        created_matrix = zeros(matrix_size);
                    
                        FD_drought_start = [];
                        FD_drought_end = [];
                    
                        % 循环处理 values 中的每个值
                        for k = 1:length(values)
                            current_value = values(k);
                    
                            drought_start = drought_starts(k);
                            drought_end = drought_ends(k);
                        
                            drought_develop1 = current_value - drought_start;
                            drought_recover1 = drought_end - current_value;
                        
                            drought_develop2 = 13 - drought_develop1;
                            drought_recover2 = 13 + drought_recover1;
                    
                            FD_drought_start = [FD_drought_start, drought_develop2];
                            FD_drought_end = [FD_drought_end, drought_recover2];
                    
                            % 检查value是否在指定范围内
                            if current_value <= 12 || current_value > 908
                                skippedCount = skippedCount + 1;
                                continue;
                            end
                        
                            % 获得前后8个索引号
                            startIndex =  current_value - 12;
                            endIndex   =  current_value + 12;
                        
                            forest_ESI = forest_ESI_data(row1, col1, startIndex:endIndex);
                    
                            forest_ESI = forest_ESI(:);
                    
                            created_matrix(1:25, k) = forest_ESI;
                        end
                    
                        created_matrix(created_matrix == 0) = nan;
                        % 计算每一行的均值
                        row_means = nanmean(created_matrix, 2);
                        FD_drought_start_mean = nanmean(FD_drought_start);
                        FD_drought_end_mean = nanmean(FD_drought_end);
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row1, col1, 1:25) = row_means;
                        forest_FD_ESI_start(row1, col1) = FD_drought_start_mean;
                        forest_FD_ESI_end(row1, col1) = FD_drought_end_mean;
                    end
                elseif isnumeric(data_variable)
                    % 如果是 double 类型的数据，则执行另一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 double 类型数据 ']);
    
                    periods_folder_path1 = fullfile(folder_path, current_folder_name);
        
                    % 获取 "summer" 文件夹下的文件列表
                    summer_periods_folder_path1 = dir(fullfile(periods_folder_path1, '*'));
                    summer_periods_folder_path2 = fullfile(periods_folder_path1, summer_periods_folder_path1(3).name);
    
                    summer_periods_files = dir(fullfile(summer_periods_folder_path2, '*'));
                    
                    % 存储文件完整路径的数组
                    summer_periods_file_paths = {};
        
                    for m = 1:numel(summer_periods_files)
                        % 检查文件名是否不是 "." 和 ".."，并且不包含数字
                        if ~strcmp(summer_periods_files(m).name, '.') && ...
                           ~strcmp(summer_periods_files(m).name, '..') && ...
                           ~any(isstrprop(summer_periods_files(m).name, 'digit'))
                       
                            % 构建完整路径并添加到数组中
                            summer_file_path = fullfile(summer_periods_folder_path2, summer_periods_files(m).name);
                            summer_periods_file_paths{end+1} = summer_file_path;
                        end
                    end
    
                    matched_files = {}; % 用于存储匹配的文件路径
    
                    % 使用 data_names{1} 构建正则表达式模式
                    pattern = ['^', regexprep(data_names{1}, '_', '\_'), '_\w+$'];
                    
                    for i = 1:numel(summer_periods_file_paths)
                        [~, filename, ~] = fileparts(summer_periods_file_paths{i}); % 获取文件名
                    
                        % 使用正则表达式模式匹配文件名
                        if ~isempty(regexp(filename, pattern, 'once'))
                            matched_files{end+1} = summer_periods_file_paths{i}; % 将匹配的文件路径添加到 matched_files 数组中
                        end
                    end
                    
                    if contains(data_names{1}, 'slow')
                        savePath1 = fullfile(output_path, 'summer', 'Slow');
                    else
                        savePath1 = fullfile(output_path, 'summer', 'Fast');
                    end
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
                    
                    periods_data1 = load(matched_files{1});
                    periods_data_names1 = fieldnames(periods_data1);
                    periods_end = periods_data1.(periods_data_names1{1});
    
                    periods_data2 = load(matched_files{2});
                    periods_data_names2 = fieldnames(periods_data2);
                    periods_start = periods_data2.(periods_data_names2{1});
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    forest_FD_ESI_start = zeros(3900, 6200);
                    forest_FD_ESI_end = zeros(3900, 6200);
    
                    % 找到非零值的位置
                    [rows2, cols2] = find(data_variable);
    
                    for i = 1:length(rows2)
                        row2 = rows2(i);
                        col2 = cols2(i);
                        value = data_variable(row2, col2);
                    
                        drought_start = periods_start(row2, col2);
                        drought_end = periods_end(row2, col2);
                    
                        drought_develop1 = value - drought_start;
                        drought_recover1 = drought_end - value;
                    
                        drought_develop2 = 13 - drought_develop1;
                        drought_recover2 = 13 + drought_recover1;
                    
                        % 检查value是否在指定范围内
                        if value <= 12 || value > 908
                            skippedCount = skippedCount + 1;
                            continue;
                        end
                    
                        % 获得前后8个索引号
                        startIndex =  value - 12;
                        endIndex =  value + 12;
                    
                        forest_ESI = forest_ESI_data(row2, col2, startIndex:endIndex);
                        % 将子集转换为列向量
                        forest_ESI = forest_ESI(:);
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row2, col2, 1:25) = forest_ESI;
                    
                        forest_FD_ESI_start(row2, col2) = drought_develop2;
                        forest_FD_ESI_end(row2, col2) = drought_recover2;
                    end         
                else
                    % 其他数据类型的处理
                    disp(['文件 ', file_paths{x}, ' 不处理']);
                end
    
                filename1 = [current_folder_name '_forest_FD_ESI_data.mat'];
                save(fullfile(savePath1, filename1), 'forest_FD_ESI_data', '-v7.3');

                filename2 = [current_folder_name '_forest_FD_ESI_start.mat'];
                save(fullfile(savePath1, filename2), 'forest_FD_ESI_start', '-v7.3');

                filename3 = [current_folder_name '_forest_FD_ESI_end.mat'];
                save(fullfile(savePath1, filename3), 'forest_FD_ESI_end', '-v7.3');
                
                fprintf('已输出 %s : forest_FD_ESI\n', current_folder_name);
            end
        elseif strcmp(current_folder_name, 'autumn')
            % 在这里执行 "autumn" 文件夹的处理操作
            files = dir(fullfile(current_folder_path, '*'));
        
            % 存储文件完整路径的数组
            file_paths = {};
            
            for m = 1:numel(files)
                % 检查文件名是否不是 "." 和 ".."，且不含有数字，并且是文件而不是文件夹
                if ~strcmp(files(m).name, '.') && ~strcmp(files(m).name, '..') && ...
                        isempty(regexp(files(m).name, '\d', 'once')) && ...
                        ~files(m).isdir
                    % 构建完整路径并存储到数组
                    file_paths{end+1} = fullfile(current_folder_path, files(m).name);
                end
            end
    
            for x = 1:numel(file_paths)
                % 加载文件
                loaded_data = load(file_paths{x});
                data_names = fieldnames(loaded_data);
                data_variable = loaded_data.(data_names{1});
                
                % 判断加载的数据类型
                if iscell(data_variable)
                   % 如果是cell数组类型的数据，则执行一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 cell 类型数据 ']);
    
                    periods_folder_path1 = fullfile(folder_path, current_folder_name);
        
                    % 获取 "autumn" 文件夹下的文件列表
                    autumn_periods_folder_path1 = dir(fullfile(periods_folder_path1, '*'));
                    autumn_periods_folder_path2 = fullfile(periods_folder_path1, autumn_periods_folder_path1(3).name);
    
                    autumn_periods_files = dir(fullfile(autumn_periods_folder_path2, '*'));
    
                    
                    % 存储文件完整路径的数组
                    autumn_periods_file_paths = {};
        
                    for m = 1:numel(autumn_periods_files)
                        % 检查文件名是否不是 "." 和 ".."，并且不包含数字
                        if ~strcmp(autumn_periods_files(m).name, '.') && ...
                           ~strcmp(autumn_periods_files(m).name, '..') && ...
                           ~any(isstrprop(autumn_periods_files(m).name, 'digit'))
                       
                            % 构建完整路径并添加到数组中
                            autumn_file_path = fullfile(autumn_periods_folder_path2, autumn_periods_files(m).name);
                            autumn_periods_file_paths{end+1} = autumn_file_path;
                        end
                    end
    
                    matched_files = {}; % 用于存储匹配的文件路径
                    % 使用 data_names{1} 构建正则表达式模式
                    pattern = ['^', regexprep(data_names{1}, '_', '\_'), '_\w+$'];
                    
                    for i = 1:numel(autumn_periods_file_paths)
                        [~, filename, ~] = fileparts(autumn_periods_file_paths{i}); % 获取文件名
                    
                        % 使用正则表达式模式匹配文件名
                        if ~isempty(regexp(filename, pattern, 'once'))
                            matched_files{end+1} = autumn_periods_file_paths{i}; % 将匹配的文件路径添加到 matched_files 数组中
                        end
                    end
    
                    savePath1 = fullfile(output_path, 'autumn', 'Mean');
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
     
                    periods_data1 = load(matched_files{1});
                    periods_data_names1 = fieldnames(periods_data1);
                    periods_end = periods_data1.(periods_data_names1{1});
    
                    periods_data2 = load(matched_files{2});
                    periods_data_names2 = fieldnames(periods_data2);
                    periods_start = periods_data2.(periods_data_names2{1});
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    forest_FD_ESI_start = zeros(3900, 6200);
                    forest_FD_ESI_end = zeros(3900, 6200);
                    
                    % 创建一个逻辑索引，表示每个单元格是否非空
                    nonEmptyCells = ~cellfun('isempty', data_variable);
                    
                    % 使用 find 函数找到非空元素的位置
                    [rows1, cols1] = find(nonEmptyCells);
                    
                    skippedCount = 0;
                    
                    for i = 1:length(rows1)
                        row1 = rows1(i);
                        col1 = cols1(i);
                        values = data_variable{row1, col1};
                        drought_starts = periods_start{row1, col1};
                        drought_ends = periods_end{row1, col1};
                    
                        % 获取 values 中的值的数量
                        num_values = numel(values);
                        
                        % 创建 25*values数量的矩阵
                        matrix_size = [25, num_values];
                        created_matrix = zeros(matrix_size);
                    
                        FD_drought_start = [];
                        FD_drought_end = [];
                    
                        % 循环处理 values 中的每个值
                        for k = 1:length(values)
                            current_value = values(k);
                    
                            drought_start = drought_starts(k);
                            drought_end = drought_ends(k);
                        
                            drought_develop1 = current_value - drought_start;
                            drought_recover1 = drought_end - current_value;
                        
                            drought_develop2 = 13 - drought_develop1;
                            drought_recover2 = 13 + drought_recover1;
                    
                            FD_drought_start = [FD_drought_start, drought_develop2];
                            FD_drought_end = [FD_drought_end, drought_recover2];
                    
                            % 检查value是否在指定范围内
                            if current_value <= 12 || current_value > 908
                                skippedCount = skippedCount + 1;
                                continue;
                            end
                        
                            % 获得前后8个索引号
                            startIndex =  current_value - 12;
                            endIndex   =  current_value + 12;
                        
                            forest_ESI = forest_ESI_data(row1, col1, startIndex:endIndex);
                    
                            forest_ESI = forest_ESI(:);
                    
                            created_matrix(1:25, k) = forest_ESI;
                        end
                    
                        created_matrix(created_matrix == 0) = nan;
                        % 计算每一行的均值
                        row_means = nanmean(created_matrix, 2);
                        FD_drought_start_mean = nanmean(FD_drought_start);
                        FD_drought_end_mean = nanmean(FD_drought_end);
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row1, col1, 1:25) = row_means;
                        forest_FD_ESI_start(row1, col1) = FD_drought_start_mean;
                        forest_FD_ESI_end(row1, col1) = FD_drought_end_mean;
                    end
                elseif isnumeric(data_variable)
                    % 如果是 double 类型的数据，则执行另一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 double 类型数据 ']);
    
                    periods_folder_path1 = fullfile(folder_path, current_folder_name);
        
                    % 获取 "autumn" 文件夹下的文件列表
                    autumn_periods_folder_path1 = dir(fullfile(periods_folder_path1, '*'));
                    autumn_periods_folder_path2 = fullfile(periods_folder_path1, autumn_periods_folder_path1(3).name);
    
                    autumn_periods_files = dir(fullfile(autumn_periods_folder_path2, '*'));
                    
                    % 存储文件完整路径的数组
                    autumn_periods_file_paths = {};
        
                    for m = 1:numel(autumn_periods_files)
                        % 检查文件名是否不是 "." 和 ".."，并且不包含数字
                        if ~strcmp(autumn_periods_files(m).name, '.') && ...
                           ~strcmp(autumn_periods_files(m).name, '..') && ...
                           ~any(isstrprop(autumn_periods_files(m).name, 'digit'))
                       
                            % 构建完整路径并添加到数组中
                            autumn_file_path = fullfile(autumn_periods_folder_path2, autumn_periods_files(m).name);
                            autumn_periods_file_paths{end+1} = autumn_file_path;
                        end
                    end
    
                    matched_files = {}; % 用于存储匹配的文件路径
    
                    % 使用 data_names{1} 构建正则表达式模式
                    pattern = ['^', regexprep(data_names{1}, '_', '\_'), '_\w+$'];
                    
                    for i = 1:numel(autumn_periods_file_paths)
                        [~, filename, ~] = fileparts(autumn_periods_file_paths{i}); % 获取文件名
                    
                        % 使用正则表达式模式匹配文件名
                        if ~isempty(regexp(filename, pattern, 'once'))
                            matched_files{end+1} = autumn_periods_file_paths{i}; % 将匹配的文件路径添加到 matched_files 数组中
                        end
                    end
                    
                    if contains(data_names{1}, 'slow')
                        savePath1 = fullfile(output_path, 'autumn', 'Slow');
                    else
                        savePath1 = fullfile(output_path, 'autumn', 'Fast');
                    end
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
                    
                    periods_data1 = load(matched_files{1});
                    periods_data_names1 = fieldnames(periods_data1);
                    periods_end = periods_data1.(periods_data_names1{1});
    
                    periods_data2 = load(matched_files{2});
                    periods_data_names2 = fieldnames(periods_data2);
                    periods_start = periods_data2.(periods_data_names2{1});
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    forest_FD_ESI_start = zeros(3900, 6200);
                    forest_FD_ESI_end = zeros(3900, 6200);
    
                    % 找到非零值的位置
                    [rows2, cols2] = find(data_variable);
    
                    for i = 1:length(rows2)
                        row2 = rows2(i);
                        col2 = cols2(i);
                        value = data_variable(row2, col2);
                    
                        drought_start = periods_start(row2, col2);
                        drought_end = periods_end(row2, col2);
                    
                        drought_develop1 = value - drought_start;
                        drought_recover1 = drought_end - value;
                    
                        drought_develop2 = 13 - drought_develop1;
                        drought_recover2 = 13 + drought_recover1;
                    
                        % 检查value是否在指定范围内
                        if value <= 12 || value > 908
                            skippedCount = skippedCount + 1;
                            continue;
                        end
                    
                        % 获得前后8个索引号
                        startIndex =  value - 12;
                        endIndex =  value + 12;
                    
                        forest_ESI = forest_ESI_data(row2, col2, startIndex:endIndex);
                        % 将子集转换为列向量
                        forest_ESI = forest_ESI(:);
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row2, col2, 1:25) = forest_ESI;
                    
                        forest_FD_ESI_start(row2, col2) = drought_develop2;
                        forest_FD_ESI_end(row2, col2) = drought_recover2;
                    end         
                else
                    % 其他数据类型的处理
                    disp(['文件 ', file_paths{x}, ' 不处理']);
                end
    
                filename1 = [current_folder_name '_forest_FD_ESI_data.mat'];
                save(fullfile(savePath1, filename1), 'forest_FD_ESI_data', '-v7.3');

                filename2 = [current_folder_name '_forest_FD_ESI_start.mat'];
                save(fullfile(savePath1, filename2), 'forest_FD_ESI_start', '-v7.3');

                filename3 = [current_folder_name '_forest_FD_ESI_end.mat'];
                save(fullfile(savePath1, filename3), 'forest_FD_ESI_end', '-v7.3');
                
                fprintf('已输出 %s : forest_FD_ESI\n', current_folder_name);
            end
        elseif strcmp(current_folder_name, 'winter')
                    % 在这里执行 "winter" 文件夹的处理操作
            files = dir(fullfile(current_folder_path, '*'));
        
            % 存储文件完整路径的数组
            file_paths = {};
            
            for m = 1:numel(files)
                % 检查文件名是否不是 "." 和 ".."，且不含有数字，并且是文件而不是文件夹
                if ~strcmp(files(m).name, '.') && ~strcmp(files(m).name, '..') && ...
                        isempty(regexp(files(m).name, '\d', 'once')) && ...
                        ~files(m).isdir
                    % 构建完整路径并存储到数组
                    file_paths{end+1} = fullfile(current_folder_path, files(m).name);
                end
            end
    
            for x = 1:numel(file_paths)
                % 加载文件
                loaded_data = load(file_paths{x});
                data_names = fieldnames(loaded_data);
                data_variable = loaded_data.(data_names{1});
                
                % 判断加载的数据类型
                if iscell(data_variable)
                   % 如果是cell数组类型的数据，则执行一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 cell 类型数据 ']);
    
                    periods_folder_path1 = fullfile(folder_path, current_folder_name);
        
                    % 获取 "winter" 文件夹下的文件列表
                    winter_periods_folder_path1 = dir(fullfile(periods_folder_path1, '*'));
                    winter_periods_folder_path2 = fullfile(periods_folder_path1, winter_periods_folder_path1(3).name);
    
                    winter_periods_files = dir(fullfile(winter_periods_folder_path2, '*'));
    
                    
                    % 存储文件完整路径的数组
                    winter_periods_file_paths = {};
        
                    for m = 1:numel(winter_periods_files)
                        % 检查文件名是否不是 "." 和 ".."，并且不包含数字
                        if ~strcmp(winter_periods_files(m).name, '.') && ...
                           ~strcmp(winter_periods_files(m).name, '..') && ...
                           ~any(isstrprop(winter_periods_files(m).name, 'digit'))
                       
                            % 构建完整路径并添加到数组中
                            winter_file_path = fullfile(winter_periods_folder_path2, winter_periods_files(m).name);
                            winter_periods_file_paths{end+1} = winter_file_path;
                        end
                    end
    
                    matched_files = {}; % 用于存储匹配的文件路径
                    % 使用 data_names{1} 构建正则表达式模式
                    pattern = ['^', regexprep(data_names{1}, '_', '\_'), '_\w+$'];
                    
                    for i = 1:numel(winter_periods_file_paths)
                        [~, filename, ~] = fileparts(winter_periods_file_paths{i}); % 获取文件名
                    
                        % 使用正则表达式模式匹配文件名
                        if ~isempty(regexp(filename, pattern, 'once'))
                            matched_files{end+1} = winter_periods_file_paths{i}; % 将匹配的文件路径添加到 matched_files 数组中
                        end
                    end
    
                    savePath1 = fullfile(output_path, 'winter', 'Mean');
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
     
                    periods_data1 = load(matched_files{1});
                    periods_data_names1 = fieldnames(periods_data1);
                    periods_end = periods_data1.(periods_data_names1{1});
    
                    periods_data2 = load(matched_files{2});
                    periods_data_names2 = fieldnames(periods_data2);
                    periods_start = periods_data2.(periods_data_names2{1});
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    forest_FD_ESI_start = zeros(3900, 6200);
                    forest_FD_ESI_end = zeros(3900, 6200);
                    
                    % 创建一个逻辑索引，表示每个单元格是否非空
                    nonEmptyCells = ~cellfun('isempty', data_variable);
                    
                    % 使用 find 函数找到非空元素的位置
                    [rows1, cols1] = find(nonEmptyCells);
                    
                    skippedCount = 0;
                    
                    for i = 1:length(rows1)
                        row1 = rows1(i);
                        col1 = cols1(i);
                        values = data_variable{row1, col1};
                        drought_starts = periods_start{row1, col1};
                        drought_ends = periods_end{row1, col1};
                    
                        % 获取 values 中的值的数量
                        num_values = numel(values);
                        
                        % 创建 25*values数量的矩阵
                        matrix_size = [25, num_values];
                        created_matrix = zeros(matrix_size);
                    
                        FD_drought_start = [];
                        FD_drought_end = [];
                    
                        % 循环处理 values 中的每个值
                        for k = 1:length(values)
                            current_value = values(k);
                    
                            drought_start = drought_starts(k);
                            drought_end = drought_ends(k);
                        
                            drought_develop1 = current_value - drought_start;
                            drought_recover1 = drought_end - current_value;
                        
                            drought_develop2 = 13 - drought_develop1;
                            drought_recover2 = 13 + drought_recover1;
                    
                            FD_drought_start = [FD_drought_start, drought_develop2];
                            FD_drought_end = [FD_drought_end, drought_recover2];
                    
                            % 检查value是否在指定范围内
                            if current_value <= 12 || current_value > 908
                                skippedCount = skippedCount + 1;
                                continue;
                            end
                        
                            % 获得前后8个索引号
                            startIndex =  current_value - 12;
                            endIndex   =  current_value + 12;
                        
                            forest_ESI = forest_ESI_data(row1, col1, startIndex:endIndex);
                    
                            forest_ESI = forest_ESI(:);
                    
                            created_matrix(1:25, k) = forest_ESI;
                        end
                    
                        created_matrix(created_matrix == 0) = nan;
                        % 计算每一行的均值
                        row_means = nanmean(created_matrix, 2);
                        FD_drought_start_mean = nanmean(FD_drought_start);
                        FD_drought_end_mean = nanmean(FD_drought_end);
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row1, col1, 1:25) = row_means;
                        forest_FD_ESI_start(row1, col1) = FD_drought_start_mean;
                        forest_FD_ESI_end(row1, col1) = FD_drought_end_mean;
                    end
                elseif isnumeric(data_variable)
                    % 如果是 double 类型的数据，则执行另一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 double 类型数据 ']);
    
                    periods_folder_path1 = fullfile(folder_path, current_folder_name);
        
                    % 获取 "winter" 文件夹下的文件列表
                    winter_periods_folder_path1 = dir(fullfile(periods_folder_path1, '*'));
                    winter_periods_folder_path2 = fullfile(periods_folder_path1, winter_periods_folder_path1(3).name);
    
                    winter_periods_files = dir(fullfile(winter_periods_folder_path2, '*'));
                    
                    % 存储文件完整路径的数组
                    winter_periods_file_paths = {};
        
                    for m = 1:numel(winter_periods_files)
                        % 检查文件名是否不是 "." 和 ".."，并且不包含数字
                        if ~strcmp(winter_periods_files(m).name, '.') && ...
                           ~strcmp(winter_periods_files(m).name, '..') && ...
                           ~any(isstrprop(winter_periods_files(m).name, 'digit'))
                       
                            % 构建完整路径并添加到数组中
                            winter_file_path = fullfile(winter_periods_folder_path2, winter_periods_files(m).name);
                            winter_periods_file_paths{end+1} = winter_file_path;
                        end
                    end
    
                    matched_files = {}; % 用于存储匹配的文件路径
    
                    % 使用 data_names{1} 构建正则表达式模式
                    pattern = ['^', regexprep(data_names{1}, '_', '\_'), '_\w+$'];
                    
                    for i = 1:numel(winter_periods_file_paths)
                        [~, filename, ~] = fileparts(winter_periods_file_paths{i}); % 获取文件名
                    
                        % 使用正则表达式模式匹配文件名
                        if ~isempty(regexp(filename, pattern, 'once'))
                            matched_files{end+1} = winter_periods_file_paths{i}; % 将匹配的文件路径添加到 matched_files 数组中
                        end
                    end
                    
                    if contains(data_names{1}, 'slow')
                        savePath1 = fullfile(output_path, 'winter', 'Slow');
                    else
                        savePath1 = fullfile(output_path, 'winter', 'Fast');
                    end
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
                    
                    periods_data1 = load(matched_files{1});
                    periods_data_names1 = fieldnames(periods_data1);
                    periods_end = periods_data1.(periods_data_names1{1});
    
                    periods_data2 = load(matched_files{2});
                    periods_data_names2 = fieldnames(periods_data2);
                    periods_start = periods_data2.(periods_data_names2{1});
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    forest_FD_ESI_start = zeros(3900, 6200);
                    forest_FD_ESI_end = zeros(3900, 6200);
    
                    % 找到非零值的位置
                    [rows2, cols2] = find(data_variable);
    
                    for i = 1:length(rows2)
                        row2 = rows2(i);
                        col2 = cols2(i);
                        value = data_variable(row2, col2);
                    
                        drought_start = periods_start(row2, col2);
                        drought_end = periods_end(row2, col2);
                    
                        drought_develop1 = value - drought_start;
                        drought_recover1 = drought_end - value;
                    
                        drought_develop2 = 13 - drought_develop1;
                        drought_recover2 = 13 + drought_recover1;
                    
                        % 检查value是否在指定范围内
                        if value <= 12 || value > 908
                            skippedCount = skippedCount + 1;
                            continue;
                        end
                    
                        % 获得前后8个索引号
                        startIndex =  value - 12;
                        endIndex =  value + 12;
                    
                        forest_ESI = forest_ESI_data(row2, col2, startIndex:endIndex);
                        % 将子集转换为列向量
                        forest_ESI = forest_ESI(:);
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row2, col2, 1:25) = forest_ESI;
                    
                        forest_FD_ESI_start(row2, col2) = drought_develop2;
                        forest_FD_ESI_end(row2, col2) = drought_recover2;
                    end         
                else
                    % 其他数据类型的处理
                    disp(['文件 ', file_paths{x}, ' 不处理']);
                end

                filename1 = [current_folder_name '_forest_FD_ESI_data.mat'];
                save(fullfile(savePath1, filename1), 'forest_FD_ESI_data', '-v7.3');

                filename2 = [current_folder_name '_forest_FD_ESI_start.mat'];
                save(fullfile(savePath1, filename2), 'forest_FD_ESI_start', '-v7.3');

                filename3 = [current_folder_name '_forest_FD_ESI_end.mat'];
                save(fullfile(savePath1, filename3), 'forest_FD_ESI_end', '-v7.3');
                
                fprintf('已输出 %s : forest_FD_ESI\n', current_folder_name);
            end
        end
    end
    fprintf('已输出 %s \n', output_path);
end