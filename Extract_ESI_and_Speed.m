clc;
clear;

% 读取栅格数据
[A,R]=readgeoraster('F:\GZW\Output\Replace.tif');

% 提供文件夹路径
folders = {
    'E:\GZW\Drought_response\Season\FD\1CI_new_1128\autumn\Mean'; 
    'E:\GZW\Drought_response\Season\FD\1CI_new_1128\spring\Mean';
    'E:\GZW\Drought_response\Season\FD\1CI_new_1128\summer\Mean';
    'E:\GZW\Drought_response\Season\FD\1CI_new_1128\winter\Mean'
};

savePath1 = 'M:\GZW\Drought_quantify\ESI&Speed\Drought_response_0417';

% 遍历每个文件夹
for i = 1:length(folders)
    folderPath = folders{i};

    [~, season_name] = fileparts(fileparts(folderPath));
    
    % 获取当前文件夹中的所有文件（包括子文件夹）
    files = dir(folderPath);

    % 存储文件路径的数组
    filePaths = {};
    
    % 遍历文件夹中的每个条目
    for j = 1:length(files)
        % 忽略 . 和 .. 以及子文件夹
        if ~files(j).isdir
            % 获取完整文件路径
            fullFilePath = fullfile(folderPath, files(j).name);
            % 将文件路径添加到 filePaths 数组中
            filePaths{end+1} = fullFilePath;
        end
    end

    % 存储文件路径的数组
    filePaths = filePaths';

%     load(filePaths{3});
% 
%     serious_ESI = forest_FD_ESI_data(:,:,13);
% 
%     output_file = sprintf('serious_ESI_%s',season_name);
% 
%     new_output_path = fullfile(savePath1, output_file);
%     geotiffwrite(new_output_path, serious_ESI, R);


    load(filePaths{1});

    output_file = sprintf('speed_%s',season_name);

    new_output_path = fullfile(savePath1, output_file);
    geotiffwrite(new_output_path, forest_FD_DESI_2D, R);
        
end