
% % 定义输入文件夹路径列表
input_folders = {
    'L:\GZW\Drought_quantify\Season_ALL_no_drought_1_month\FD',...
    'L:\GZW\Drought_quantify\Season_N&P_no_drought_1_month\Natural\FD',...
    'L:\GZW\Drought_quantify\Season_N&P_no_drought_1_month\Planted\FD'
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
    
                    savePath1 = fullfile(output_path, 'spring');
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    
                    % 创建一个逻辑索引，表示每个单元格是否非空
                    nonEmptyCells = ~cellfun('isempty', data_variable);
                    
                    % 使用 find 函数找到非空元素的位置
                    [rows1, cols1] = find(nonEmptyCells);
                    
                    skippedCount = 0;
                    
                    for i = 1:length(rows1)
                        row1 = rows1(i);
                        col1 = cols1(i);
                        values = data_variable{row1, col1};
                    
                        % 获取 values 中的值的数量
                        num_values = numel(values);
                        
                        % 创建 25*values数量的矩阵
                        matrix_size = [25, num_values];
                        created_matrix = zeros(matrix_size);
                    
                        % 循环处理 values 中的每个值
                        for k = 1:length(values)
                            current_value = values(k);
                    
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
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row1, col1, 1:25) = row_means;
                    end
                end

                filename1 = [current_folder_name '_forest_FD_ESI_data.mat'];
                save(fullfile(savePath1, filename1), 'forest_FD_ESI_data', '-v7.3');
                
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
    
                    savePath1 = fullfile(output_path, 'summer');
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    
                    % 创建一个逻辑索引，表示每个单元格是否非空
                    nonEmptyCells = ~cellfun('isempty', data_variable);
                    
                    % 使用 find 函数找到非空元素的位置
                    [rows1, cols1] = find(nonEmptyCells);
                    
                    skippedCount = 0;
                    
                    for i = 1:length(rows1)
                        row1 = rows1(i);
                        col1 = cols1(i);
                        values = data_variable{row1, col1};
                    
                        % 获取 values 中的值的数量
                        num_values = numel(values);
                        
                        % 创建 25*values数量的矩阵
                        matrix_size = [25, num_values];
                        created_matrix = zeros(matrix_size);
                    
                        % 循环处理 values 中的每个值
                        for k = 1:length(values)
                            current_value = values(k);
                    
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
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row1, col1, 1:25) = row_means;
                    end
                end
                
                filename1 = [current_folder_name '_forest_FD_ESI_data.mat'];
                save(fullfile(savePath1, filename1), 'forest_FD_ESI_data', '-v7.3');
                
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
    
                    savePath1 = fullfile(output_path, 'autumn');
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    
                    % 创建一个逻辑索引，表示每个单元格是否非空
                    nonEmptyCells = ~cellfun('isempty', data_variable);
                    
                    % 使用 find 函数找到非空元素的位置
                    [rows1, cols1] = find(nonEmptyCells);
                    
                    skippedCount = 0;
                    
                    for i = 1:length(rows1)
                        row1 = rows1(i);
                        col1 = cols1(i);
                        values = data_variable{row1, col1};
                    
                        % 获取 values 中的值的数量
                        num_values = numel(values);
                        
                        % 创建 25*values数量的矩阵
                        matrix_size = [25, num_values];
                        created_matrix = zeros(matrix_size);
                    
                        % 循环处理 values 中的每个值
                        for k = 1:length(values)
                            current_value = values(k);
                    
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
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row1, col1, 1:25) = row_means;
                    end
                end

                filename1 = [current_folder_name '_forest_FD_ESI_data.mat'];
                save(fullfile(savePath1, filename1), 'forest_FD_ESI_data', '-v7.3');
                
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
    
                    savePath1 = fullfile(output_path, 'winter');
    
                    % 检查保存路径是否存在，不存在则创建
                    if ~isfolder(savePath1)
                        mkdir(savePath1); % 创建保存路径
                    end
    
                    forest_FD_ESI_data = zeros(3900, 6200, 25);%   ????
                    
                    % 创建一个逻辑索引，表示每个单元格是否非空
                    nonEmptyCells = ~cellfun('isempty', data_variable);
                    
                    % 使用 find 函数找到非空元素的位置
                    [rows1, cols1] = find(nonEmptyCells);
                    
                    skippedCount = 0;
                    
                    for i = 1:length(rows1)
                        row1 = rows1(i);
                        col1 = cols1(i);
                        values = data_variable{row1, col1};
                    
                        % 获取 values 中的值的数量
                        num_values = numel(values);
                        
                        % 创建 25*values数量的矩阵
                        matrix_size = [25, num_values];
                        created_matrix = zeros(matrix_size);
                    
                        % 循环处理 values 中的每个值
                        for k = 1:length(values)
                            current_value = values(k);
                    
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
                    
                        % 将 forest_ESI 放入对应索引的位置
                        forest_FD_ESI_data(row1, col1, 1:25) = row_means;
                    end
                end
    
                filename1 = [current_folder_name '_forest_FD_ESI_data.mat'];
                save(fullfile(savePath1, filename1), 'forest_FD_ESI_data', '-v7.3');
                
                fprintf('已输出 %s : forest_FD_ESI\n', current_folder_name);
            end
        end
    end
    fprintf('已输出 %s \n', output_path);
end