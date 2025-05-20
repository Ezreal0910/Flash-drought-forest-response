clc;
clear;

% 读取栅格数据
[A,R]=readgeoraster('F:\GZW\Output\Replace.tif');

% 定义输入文件夹路径列表
input_folders = {
    'L:\GZW\Drought_quantify\Season_ALL\FD', ...
    'L:\GZW\Drought_quantify\Season_ALL_1_month\FD',...
    'L:\GZW\Drought_quantify\Season_N&P\Natural\FD',...
    'L:\GZW\Drought_quantify\Season_N&P\Planted\FD',...
    'L:\GZW\Drought_quantify\Season_N&P_1_month\Natural\FD',...
    'L:\GZW\Drought_quantify\Season_N&P_1_month\Planted\FD'
};

% 定义新的输出子文件夹名称
target_parent = 'E:\GZW\Drought_response\';

% 新的子文件夹名称
new_subfolder = 'Drought_distributed2';

% 循环处理每个输入文件夹
for f = 1:length(input_folders)
    folder_path = input_folders{f};

    % 判断文件夹路径中是否包含 "all"
    if contains(folder_path, 'all')
        % 获取文件夹中的所有文件
        files = dir(fullfile(folder_path, '*'));
        files = files(~[files.isdir]);  % 仅保留文件

        % 存储文件完整路径的数组
        file_paths = {};

        for m = 1:numel(files)
            % 检查文件名是否不是 "." 和 ".."，且文件名不包含 'FD' 和 'SD'
            if ~strcmp(files(m).name, '.') && ~strcmp(files(m).name, '..') && ...
                    isempty(regexp(files(m).name, '\d.*\d', 'once')) && ...
                    ~contains(files(m).name, 'FD') && ...
                    ~contains(files(m).name, 'SD') && ...
                    ~files(m).isdir
                % 构建完整路径并存储到数组
                file_paths{end+1} = fullfile(folder_path, files(m).name);
            end
        end

        for x = 1:numel(file_paths)
            file_path = file_paths{x};
        
            % 获取文件名
            [~, file_name, ext] = fileparts(file_path);

            % 初始化新文件名
            tree_name = '';

            % 检查文件名是否包含数字 0 或 1
            if contains(file_name, '0')
                tree_name = 'Natural';
            elseif contains(file_name, '1')
                tree_name = 'Planted';
            end

            full_file_name = strcat(file_name, ext);

            % 加载文件
            loaded_data = load(file_paths{x});
            data_names = fieldnames(loaded_data);
            data_variable = loaded_data.(data_names{1});
            
            % 判断加载的数据类型
            if iscell(data_variable)
               % 如果是cell数组类型的数据，则执行一种处理
                disp(['文件 ', file_paths{x}, ' 包含 cell 类型数据 ']);
                
                % 使用 cellfun 函数计算每个单元格中的元素数量
                drought_distributed = cellfun(@numel, data_variable);

                drought_distributed(drought_distributed == 0) = 99;
                
                % 检查文件名并构建新的输出文件名
                if contains(full_file_name, 'Flash_Drought')
                    if ~isempty(tree_name)
                        new_file_name = strcat('ALL_', tree_name,'_FD', '.tif');
                    else
                        new_file_name = 'ALL_FD.tif';
                    end
                elseif contains(full_file_name, 'Slow_Drought')
                    if ~isempty(tree_name)
                        new_file_name = strcat('ALL_', tree_name,'_SD', '.tif');
                    else
                        new_file_name = 'ALL_SD.tif';
                    end
                end
                
                % 构建新的输出路径
                new_output_path = fullfile(target_parent, new_subfolder, new_file_name);
                
                % 写回栅格数据
                geotiffwrite(new_output_path, drought_distributed, R);
                fprintf('已输出 %s \n', new_output_path);
            end
        end
    else
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

        % 获取相对于目标父路径的子文件夹名
        relative_path = extractAfter(folder_path, target_parent);
        relative_path = strrep(relative_path, '\', '_'); % 将 '\' 替换为 '_'

        for j = 1:length(folder_names)
            current_folder_name = folder_names{j};
    
            % 检查是否为 '1CI'，如果是则跳过当前循环迭代
            if strcmp(current_folder_name, '1CI')
                continue;  % 跳过当前迭代
            end
    
            current_folder_path = fullfile(folder_path, current_folder_name);
            
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
                file_path = file_paths{x};
        
                % 获取文件名
                [~, file_name, ext] = fileparts(file_path);
                full_file_name = strcat(file_name, ext);
            
                % 检查文件名中是否包含 'FD' 或 'SD'
                if contains(full_file_name, 'FD') || contains(full_file_name, 'SD')
                    continue; % 跳过该文件
                end

                % 加载文件
                loaded_data = load(file_paths{x});
                data_names = fieldnames(loaded_data);
                data_variable = loaded_data.(data_names{1});
                
                % 判断加载的数据类型
                if iscell(data_variable)
                   % 如果是cell数组类型的数据，则执行一种处理
                    disp(['文件 ', file_paths{x}, ' 包含 cell 类型数据 ']);
                    
                    % 使用 cellfun 函数计算每个单元格中的元素数量
                    drought_distributed = cellfun(@numel, data_variable);

                    drought_distributed(drought_distributed == 0) = 99;
                    
                    % 构建新的输出文件名
                    new_file_name = [relative_path, '_', current_folder_name, '.tif'];
                    
                    % 构建新的输出路径
                    new_output_path = fullfile(target_parent, new_subfolder, new_file_name);
                    
                    % 写回栅格数据
                    geotiffwrite(new_output_path, drought_distributed, R);
                    fprintf('已输出 %s \n', new_output_path);
                end
            end
        end
    end
end
