clc;
clear;

Huanan_FD_495 = "E:\GZW\Drought_response\Region\Huanan\Dinghushan_region\Region_Huanan_FD_35_495.tif";
Huanan_FD_495_location = imread(Huanan_FD_495);

% 定义区域范围
region_rows_Huanan = 2401:3400;
region_cols_Huanan = 3201:4300;

% 提取位置数据
Huanan_FD_495_location_extracted = Huanan_FD_495_location(region_rows_Huanan, region_cols_Huanan);

% 找到有值的位置
non_zero_indices_Huanan_495 = find(Huanan_FD_495_location_extracted);

% 创建一个新的位置数据矩阵
combined_location_Huanan = zeros(size(Huanan_FD_495_location_extracted));

combined_location_Huanan(non_zero_indices_Huanan_495) = Huanan_FD_495_location_extracted(non_zero_indices_Huanan_495);

% 要计算的置信水平数组
confidence_levels = 5:10:95;

% 转换为累积概率
cumulative_probabilities = (100 + confidence_levels) / 200;

% 计算对应的 z 值
z_values = norminv(cumulative_probabilities);

% 加载 Huanan 数据
load("G:\GZW\Forest_data\Huanan\forest_FD_Events_Huanan.mat")

load("G:\GZW\Forest_data\Huanan\forest_sif_Huanan.mat")
load("G:\GZW\Forest_data\Huanan\forest_LAI_Huanan.mat")

% 找到非零值的位置
[rows_Huanan, cols_Huanan] = find(combined_location_Huanan);
% 找到非零值的位置
[rows_Huanan2, cols_Huanan2] = find(Huanan_FD_495_location_extracted);

% 创建另一个矩阵，用作折线图对比
forest_FD_sif_data_Huanan = zeros(size(forest_FD_Events_Huanan));
forest_FD_LAI_data_Huanan = zeros(size(forest_FD_Events_Huanan));

extract_event = 495;

% 获得前后12个索引号
startIndex =  extract_event - 12;
endIndex   =  extract_event + 12;

% 遍历非零元素的行和列索引
for n = 1:length(rows_Huanan)
    row_Huanan = rows_Huanan(n);
    col_Huanan = cols_Huanan(n);

    point_sif_3D_Huanan = forest_sif_Huanan(row_Huanan, col_Huanan, : );
    point_sif_1D_Huanan = point_sif_3D_Huanan(:);

    point_LAI_3D_Huanan = forest_LAI_Huanan(row_Huanan, col_Huanan, : );
    point_LAI_1D_Huanan = point_LAI_3D_Huanan(:);

    forest_sif_Huanan_drought = point_sif_1D_Huanan(startIndex:endIndex);
    forest_FD_sif_data_Huanan(row_Huanan, col_Huanan, 1:25) = forest_sif_Huanan_drought;

    forest_LAI_Huanan_drought = point_LAI_1D_Huanan(startIndex:endIndex);
    forest_FD_LAI_data_Huanan(row_Huanan, col_Huanan, 1:25) = forest_LAI_Huanan_drought;
end

forest_FD_sif_data_Huanan(forest_FD_sif_data_Huanan == 0) = nan;
forest_FD_LAI_data_Huanan(forest_FD_LAI_data_Huanan == 0) = nan;

% 预分配数组用于存储结果
forest_FD_sif_nanmean_1D = zeros(25, 1);
forest_FD_sif_nanstd_1D = zeros(25, 1);

forest_FD_sif_standard_error_1D = zeros(25, 1);
forest_FD_sif_count = zeros(25, 1);

% 预分配数组用于存储结果
forest_FD_sif_confidence_intervals_upper = zeros(25, length(confidence_levels));
forest_FD_sif_confidence_intervals_lower = zeros(25, length(confidence_levels));

for slice_index = 1:25
    % 逐个提取切片
    current_slice = forest_FD_sif_data_Huanan(:, :, slice_index);

    % 在切片上计算均值
    forest_FD_sif_nanmean_1D(slice_index) = nanmean(current_slice, 'all');

    % 在切片上计算标准差
    forest_FD_sif_nanstd_1D(slice_index) = nanstd(current_slice, 0, 'all');

    % 在切片上计算标准误差
    num_observation = sum(~isnan(current_slice(:)));  % 计算非 NaN 值的观测数量
    forest_FD_sif_count(slice_index) = num_observation;

    forest_FD_sif_standard_error_1D(slice_index) = forest_FD_sif_nanstd_1D(slice_index) / sqrt(num_observation);

    % 计算不同置信度下的置信区间
    for level_index = 1:length(confidence_levels)
        % 计算置信区间的上限和下限
        upper_limit = forest_FD_sif_nanmean_1D(slice_index) + z_values(level_index) * forest_FD_sif_standard_error_1D(slice_index);
        lower_limit = forest_FD_sif_nanmean_1D(slice_index) - z_values(level_index) * forest_FD_sif_standard_error_1D(slice_index);
        
        forest_FD_sif_confidence_intervals_upper(slice_index, level_index) = upper_limit;
        forest_FD_sif_confidence_intervals_lower(slice_index, level_index) = lower_limit;
    end
end

output_path = "E:\GZW\Drought_response\Region\Contrast\FD_495_4";

filename1 = 'FD_495_forest_sif_mean.csv';
writematrix(forest_FD_sif_nanmean_1D, fullfile(output_path, filename1));

filename2 = 'FD_495_forest_sif_ci_upper.csv';
writematrix(forest_FD_sif_confidence_intervals_upper, fullfile(output_path, filename2));

filename3 = 'FD_495_forest_sif_ci_lower.csv';
writematrix(forest_FD_sif_confidence_intervals_lower, fullfile(output_path, filename3));

filename4 = 'FD_495_forest_sif_standard_error.csv';
writematrix(forest_FD_sif_standard_error_1D, fullfile(output_path, filename4));

filename5 = 'FD_495_forest_sif_std.csv';
writematrix(forest_FD_sif_nanstd_1D, fullfile(output_path, filename5));

filename6 = 'FD_495_forest_sif_count.csv';
writematrix(forest_FD_sif_count, fullfile(output_path, filename6));

% 预分配数组用于存储结果
forest_FD_LAI_nanmean_1D = zeros(25, 1);
forest_FD_LAI_nanstd_1D = zeros(25, 1);

forest_FD_LAI_standard_error_1D = zeros(25, 1);
forest_FD_LAI_count = zeros(25, 1);

% 预分配数组用于存储结果
forest_FD_LAI_confidence_intervals_upper = zeros(25, length(confidence_levels));
forest_FD_LAI_confidence_intervals_lower = zeros(25, length(confidence_levels));

for slice_index = 1:25
    % 逐个提取切片
    current_slice = forest_FD_LAI_data_Huanan(:, :, slice_index);

    % 在切片上计算均值
    forest_FD_LAI_nanmean_1D(slice_index) = nanmean(current_slice, 'all');

    % 在切片上计算标准差
    forest_FD_LAI_nanstd_1D(slice_index) = nanstd(current_slice, 0, 'all');

    % 在切片上计算标准误差
    num_observation = sum(~isnan(current_slice(:)));  % 计算非 NaN 值的观测数量
    forest_FD_LAI_count(slice_index) = num_observation;

    forest_FD_LAI_standard_error_1D(slice_index) = forest_FD_LAI_nanstd_1D(slice_index) / sqrt(num_observation);

    % 计算不同置信度下的置信区间
    for level_index = 1:length(confidence_levels)
        % 计算置信区间的上限和下限
        upper_limit = forest_FD_LAI_nanmean_1D(slice_index) + z_values(level_index) * forest_FD_LAI_standard_error_1D(slice_index);
        lower_limit = forest_FD_LAI_nanmean_1D(slice_index) - z_values(level_index) * forest_FD_LAI_standard_error_1D(slice_index);
        
        forest_FD_LAI_confidence_intervals_upper(slice_index, level_index) = upper_limit;
        forest_FD_LAI_confidence_intervals_lower(slice_index, level_index) = lower_limit;
    end
end

output_path = "E:\GZW\Drought_response\Region\Contrast\FD_495_4";

filename1 = 'FD_495_forest_LAI_mean.csv';
writematrix(forest_FD_LAI_nanmean_1D, fullfile(output_path, filename1));

filename2 = 'FD_495_forest_LAI_ci_upper.csv';
writematrix(forest_FD_LAI_confidence_intervals_upper, fullfile(output_path, filename2));

filename3 = 'FD_495_forest_LAI_ci_lower.csv';
writematrix(forest_FD_LAI_confidence_intervals_lower, fullfile(output_path, filename3));

filename4 = 'FD_495_forest_LAI_standard_error.csv';
writematrix(forest_FD_LAI_standard_error_1D, fullfile(output_path, filename4));

filename5 = 'FD_495_forest_LAI_std.csv';
writematrix(forest_FD_LAI_nanstd_1D, fullfile(output_path, filename5));

filename6 = 'FD_495_forest_LAI_count.csv';
writematrix(forest_FD_LAI_count, fullfile(output_path, filename6));

disp('所有文件已保存');