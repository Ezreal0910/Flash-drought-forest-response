clc;
clear;

Huanan_FD_178 = "E:\GZW\Drought_response\Region\Huanan\Region_Huanan_FD_5_178.tif";
Huanan_FD_178_location = imread(Huanan_FD_178);

% 定义区域范围
region_rows_Huanan = 2401:3400;
region_cols_Huanan = 3201:4300;

% 提取位置数据
Huanan_FD_178_location_extracted = Huanan_FD_178_location(region_rows_Huanan, region_cols_Huanan);

% 找到有值的位置
non_zero_indices_Huanan_178 = find(Huanan_FD_178_location_extracted);

% 创建一个新的位置数据矩阵
combined_location_Huanan = zeros(size(Huanan_FD_178_location_extracted));

% 将第一个位置数据的有值位置复制到新矩阵的相应位置
combined_location_Huanan(non_zero_indices_Huanan_178) = Huanan_FD_178_location_extracted(non_zero_indices_Huanan_178);

Dongnan_FD_178 = "E:\GZW\Drought_response\Region\Dongnan\Region_Dongnan_FD_2_178.tif";
Dongnan_FD_178_location = imread(Dongnan_FD_178);

% 定义区域范围
region_rows_Dongnan = 2401:3400;
region_cols_Dongnan = 4301:4900;

% 提取位置数据
Dongnan_FD_178_location_extracted = Dongnan_FD_178_location(region_rows_Dongnan, region_cols_Dongnan);

% 将所有有值的位置的值改为1
Dongnan_FD_178_location_extracted(Dongnan_FD_178_location_extracted > 0) = 1;


% % 读取栅格数据
% [A,R]=readgeoraster('F:\GZW\Output\Replace.tif');
% 
% C = zeros(size(A));
% 
% % 将区域数据赋给相应区域
% C(region_rows_Huanan, region_cols_Huanan) = combined_location_Huanan;
% C(region_rows_Dongnan, region_cols_Dongnan) = Dongnan_FD_178_location_extracted;
% 
% % 写回栅格数据
% output_filename = "E:\GZW\Drought_response\Region\FD_178.tif";
% geotiffwrite(output_filename, C, R);

% 要计算的置信水平数组
confidence_levels = 5:10:95;

% 转换为累积概率
cumulative_probabilities = (100 + confidence_levels) / 200;

% 计算对应的 z 值
z_values = norminv(cumulative_probabilities);

% 加载 Huanan 数据
load("G:\GZW\Forest_data\Huanan\forest_FD_Events_Huanan.mat")
load("G:\GZW\Forest_data\Huanan\forest_GPP_data_Huanan.mat")

% 找到非零值的位置
[rows_Huanan, cols_Huanan] = find(combined_location_Huanan);

% 创建另一个矩阵，用作折线图对比
forest_FD_GPP_data_Huanan = zeros(size(forest_FD_Events_Huanan));

extract_event = 178;

% 获得前后12个索引号
startIndex =  extract_event - 12;
endIndex   =  extract_event + 12;

% 遍历非零元素的行和列索引
for n = 1:length(rows_Huanan)
    row_Huanan = rows_Huanan(n);
    col_Huanan = cols_Huanan(n);

    point_GPP_3D_Huanan = forest_GPP_data_Huanan(row_Huanan, col_Huanan, : );
    point_GPP_1D_Huanan = point_GPP_3D_Huanan(:);

    forest_GPP_data_Huanan_drought = point_GPP_1D_Huanan(startIndex:endIndex);
    forest_FD_GPP_data_Huanan(row_Huanan, col_Huanan, 1:25) = forest_GPP_data_Huanan_drought;
end

% 加载 Dongnan 数据
load("G:\GZW\Forest_data\Dongnan\forest_FD_Events_Dongnan.mat")
load("G:\GZW\Forest_data\Dongnan\forest_GPP_data_Dongnan.mat")

% 找到非零值的位置
[rows_Dongnan, cols_Dongnan] = find(Dongnan_FD_178_location_extracted);

% 创建另一个矩阵，用作折线图对比
forest_FD_GPP_data_Dongnan = zeros(size(forest_FD_Events_Dongnan));

% 遍历非零元素的行和列索引
for n = 1:length(rows_Dongnan)
    row_Dongnan = rows_Dongnan(n);
    col_Dongnan = cols_Dongnan(n);

    point_GPP_3D_Dongnan = forest_GPP_data_Dongnan(row_Dongnan, col_Dongnan, : );
    point_GPP_1D_Dongnan = point_GPP_3D_Dongnan(:);

    forest_GPP_data_Dongnan_drought = point_GPP_1D_Dongnan(startIndex:endIndex);
    forest_FD_GPP_data_Dongnan(row_Dongnan, col_Dongnan, 1:25) = forest_GPP_data_Dongnan_drought;
end

% 沿第2维度拼接两个矩阵
result_forest_FD_GPP_data = cat(2, forest_FD_GPP_data_Huanan, forest_FD_GPP_data_Dongnan);

result_forest_FD_GPP_data(result_forest_FD_GPP_data == 0) = nan;

% 预分配数组用于存储结果
forest_FD_GPP_nanmean_1D = zeros(25, 1);
forest_FD_GPP_nanstd_1D = zeros(25, 1);

forest_FD_GPP_standard_error_1D = zeros(25, 1);
forest_FD_GPP_count = zeros(25, 1);

% 预分配数组用于存储结果
forest_FD_GPP_confidence_intervals_upper = zeros(25, length(confidence_levels));
forest_FD_GPP_confidence_intervals_lower = zeros(25, length(confidence_levels));

for slice_index = 1:25
    % 逐个提取切片
    current_slice = result_forest_FD_GPP_data(:, :, slice_index);

    % 在切片上计算均值
    forest_FD_GPP_nanmean_1D(slice_index) = nanmean(current_slice, 'all');

    % 在切片上计算标准差
    forest_FD_GPP_nanstd_1D(slice_index) = nanstd(current_slice, 0, 'all');

    % 在切片上计算标准误差
    num_observation = sum(~isnan(current_slice(:)));  % 计算非 NaN 值的观测数量
    forest_FD_GPP_count(slice_index) = num_observation;

    forest_FD_GPP_standard_error_1D(slice_index) = forest_FD_GPP_nanstd_1D(slice_index) / sqrt(num_observation);

    % 计算不同置信度下的置信区间
    for level_index = 1:length(confidence_levels)
        % 计算置信区间的上限和下限
        upper_limit = forest_FD_GPP_nanmean_1D(slice_index) + z_values(level_index) * forest_FD_GPP_standard_error_1D(slice_index);
        lower_limit = forest_FD_GPP_nanmean_1D(slice_index) - z_values(level_index) * forest_FD_GPP_standard_error_1D(slice_index);
        
        forest_FD_GPP_confidence_intervals_upper(slice_index, level_index) = upper_limit;
        forest_FD_GPP_confidence_intervals_lower(slice_index, level_index) = lower_limit;
    end
end

output_path = "E:\GZW\Drought_response\Region\Contrast\FD_178_2";

filename1 = 'FD_178_forest_GPP_mean.csv';
writematrix(forest_FD_GPP_nanmean_1D, fullfile(output_path, filename1));

filename2 = 'FD_178_forest_GPP_ci_upper.csv';
writematrix(forest_FD_GPP_confidence_intervals_upper, fullfile(output_path, filename2));

filename3 = 'FD_178_forest_GPP_ci_lower.csv';
writematrix(forest_FD_GPP_confidence_intervals_lower, fullfile(output_path, filename3));

filename4 = 'FD_178_forest_GPP_standard_error.csv';
writematrix(forest_FD_GPP_standard_error_1D, fullfile(output_path, filename4));

filename5 = 'FD_178_forest_GPP_std.csv';
writematrix(forest_FD_GPP_nanstd_1D, fullfile(output_path, filename5));

filename6 = 'FD_178_forest_GPP_count.csv';
writematrix(forest_FD_GPP_count, fullfile(output_path, filename6));

disp('所有文件已保存');