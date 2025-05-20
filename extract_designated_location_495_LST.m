clc;
clear;

Huanan_FD_495 = "E:\GZW\Drought_response\Region\Huanan\Dinghushan_region\Region_Huanan_FD_35_495.tif";
Huanan_FD_495_location = imread(Huanan_FD_495);

% 定义区域范围
region_rows_Huanan = 2401:3400;
region_cols_Huanan = 3201:4300;

% 提取位置数据
Huanan_FD_495_location_extracted = Huanan_FD_495_location(region_rows_Huanan, region_cols_Huanan);

% 转为0-1：有值设为1
Huanan_binary = Huanan_FD_495_location_extracted > 0;

% basePath = 'N:\GZW\ERA5\LST_daily_mean_tif_reproject';
basePath = 'N:\GZW\LST_resample';
filePaths = []; % 存储最终文件路径的单元格数组

% 预定义矩阵尺寸
matrix2004 = zeros(3900, 6200, 189); % 单精度浮点类型
matrix2005 = zeros(3900, 6200, 8);

% 处理 2004 年数据
year2004_path = fullfile(basePath, '2011');
tif2004 = dir(fullfile(year2004_path, '*.tif'));
fileList2004 = {tif2004.name};
[~, idx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), fileList2004));
sorted2004 = fileList2004(idx);
selected2004 = sorted2004(177:end);

for i = 1:length(selected2004)
    imgPath = fullfile(year2004_path, selected2004{i});
    matrix2004(:,:,i) = single(imread(imgPath));
    fprintf('2004年填充进度: %d/%d (%.1f%%)\n', i, length(selected2004), i/length(selected2004)*100);
end


% 处理 2005 年数据
year2005_path = fullfile(basePath, '2012');
tif2005 = dir(fullfile(year2005_path, '*.tif'));
fileList2005 = {tif2005.name};
[~, idx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), fileList2005));
sorted2005 = fileList2005(idx);

for i = 1:8
    imgPath = fullfile(year2005_path, sorted2005{i});
    matrix2005(:,:,i) = single(imread(imgPath));
    fprintf('2005年填充进度: %d/48 (%.1f%%)\n', i, i/48 * 100);
end

% 替换无效值
matrix2004(matrix2004 == -9999) = NaN;
matrix2005(matrix2005 == -9999) = NaN;

% 按 8 个一组求均值（使用循环）
nGroups2004 = ceil(189 / 8);
matrix2004_avg = zeros(3900, 6200, nGroups2004);
for i = 1:nGroups2004
    startIdx = (i - 1) * 8 + 1;
    endIdx = min(startIdx + 7, 189);
    matrix2004_avg(:,:,i) = mean(matrix2004(:,:,startIdx:endIdx), 3, 'omitnan');
end

nGroups2005 = ceil(8 / 8);
matrix2005_avg = zeros(3900, 6200, nGroups2005);
for i = 1:nGroups2005
    startIdx = (i - 1) * 8 + 1;
    endIdx = min(startIdx + 7, 8);
    matrix2005_avg(:,:,i) = mean(matrix2005(:,:,startIdx:endIdx), 3, 'omitnan');
end

% 按顺序堆叠两个矩阵
matrix_avg_combined = cat(3, matrix2004_avg, matrix2005_avg);

clear matrix2004 matrix2005

% 提取指定区域的数据
matrix_avg_Huanan = matrix_avg_combined(2401:3400, 3201:4300, :);

% 找到非零值的位置
[rows_Huanan, cols_Huanan] = find(Huanan_binary);

% 转换为线性索引
linear_idx = sub2ind([1000, 1100], rows_Huanan, cols_Huanan);

% 提取层数信息
[~, ~, numLayers] = size(matrix_avg_Huanan);

% 初始化结果
matrix_avg_Huanan_extracted = zeros(length(linear_idx), numLayers);

% 循环提取每一层
for k = 1:numLayers
    layer = matrix_avg_Huanan(:, :, k);
    matrix_avg_Huanan_extracted(:, k) = layer(linear_idx);
end

% 要计算的置信水平数组
confidence_levels = 5:10:95;

% 转换为累积概率
cumulative_probabilities = (100 + confidence_levels) / 200;

% 计算对应的 z 值
z_values = norminv(cumulative_probabilities);

% 将 0 值替换为 NaN
matrix_avg_Huanan_extracted(matrix_avg_Huanan_extracted == 0) = NaN;

% 预分配数组用于存储结果
forest_FD_T_nanmean_1D = zeros(25, 1);
forest_FD_T_nanstd_1D = zeros(25, 1);

forest_FD_T_standard_error_1D = zeros(25, 1);
forest_FD_T_count = zeros(25, 1);

% 预分配数组用于存储结果
forest_FD_T_confidence_intervals_upper = zeros(25, length(confidence_levels));
forest_FD_T_confidence_intervals_lower = zeros(25, length(confidence_levels));

for slice_index = 1:25
    % 逐个提取切片
    current_slice = matrix_avg_Huanan_extracted(:, slice_index);

    % 在切片上计算均值
    forest_FD_T_nanmean_1D(slice_index) = nanmean(current_slice, 'all');

    % 在切片上计算标准差
    forest_FD_T_nanstd_1D(slice_index) = nanstd(current_slice, 0, 'all');

    % 在切片上计算标准误差
    num_observation = sum(~isnan(current_slice(:)));  % 计算非 NaN 值的观测数量
    forest_FD_T_count(slice_index) = num_observation;

    forest_FD_T_standard_error_1D(slice_index) = forest_FD_T_nanstd_1D(slice_index) / sqrt(num_observation);

    % 计算不同置信度下的置信区间
    for level_index = 1:length(confidence_levels)
        % 计算置信区间的上限和下限
        upper_limit = forest_FD_T_nanmean_1D(slice_index) + z_values(level_index) * forest_FD_T_standard_error_1D(slice_index);
        lower_limit = forest_FD_T_nanmean_1D(slice_index) - z_values(level_index) * forest_FD_T_standard_error_1D(slice_index);
        
        forest_FD_T_confidence_intervals_upper(slice_index, level_index) = upper_limit;
        forest_FD_T_confidence_intervals_lower(slice_index, level_index) = lower_limit;
    end
end

output_path = "E:\GZW\Drought_response\Region\Contrast\FD_495_4";

filename1 = 'FD_495_LST_mean.csv';
writematrix(forest_FD_T_nanmean_1D, fullfile(output_path, filename1));

filename2 = 'FD_495_LST_ci_upper.csv';
writematrix(forest_FD_T_confidence_intervals_upper, fullfile(output_path, filename2));

filename3 = 'FD_495_LST_ci_lower.csv';
writematrix(forest_FD_T_confidence_intervals_lower, fullfile(output_path, filename3));

filename4 = 'FD_495_LST_standard_error.csv';
writematrix(forest_FD_T_standard_error_1D, fullfile(output_path, filename4));

filename5 = 'FD_495_LST_std.csv';
writematrix(forest_FD_T_nanstd_1D, fullfile(output_path, filename5));

filename6 = 'FD_495_LST_count.csv';
writematrix(forest_FD_T_count, fullfile(output_path, filename6));

disp('所有文件已保存');
