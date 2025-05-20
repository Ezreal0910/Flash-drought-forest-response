clc;
clear;

result_filename = "N:\GZW\站点数据_1217\千烟洲\Result\2005_new\output.txt";

% 使用反斜杠分割路径
parts = strsplit(result_filename, '\');

% 提取目标部分
site_name = parts{4};

input_data_file = "N:\GZW\站点数据_1217\千烟洲\Result\2005_output_data_new.txt";

meteoro_data_file = "N:\GZW\站点数据_1217\千烟洲\气象数据\日尺度\2005年千烟洲气象日统计数据.xlsx";
meteoro_data = readtable(meteoro_data_file, 'VariableNamingRule', 'preserve');
P_data = meteoro_data.('日降水量');

uWUE_data_filename = "N:\GZW\站点数据_1217\千烟洲\E&T\output_daily_2005_new.xlsx";
uWUE_data = readtable(uWUE_data_filename, 'VariableNamingRule', 'preserve');
uWUE_T_data = uWUE_data.('uWUE_T');

year = meteoro_data(1,1);
year = table2array(year);

% 使用 readtable 读取数据，自动处理分隔符
original_data = readtable(input_data_file, 'VariableNamingRule', 'preserve');

NEE_data = original_data.('NEE');
Tair_data = original_data.('Tair');
VPD_data = original_data.('VPD');

LE_data = original_data.('LE');

% 使用 readtable 读取数据，自动处理分隔符
result_data = readtable(result_filename, 'VariableNamingRule', 'preserve');

GPP_calculate_columns = result_data(:, contains(result_data.Properties.VariableNames, 'GPP'));
GPP_DT_data = GPP_calculate_columns(:,4);
GPP_DT_data = table2array(GPP_DT_data);

GPP_f_data = GPP_calculate_columns(:,1);
GPP_f_data = table2array(GPP_f_data);

Reco_calculate_columns = result_data(:, contains(result_data.Properties.VariableNames, 'Reco'));
Reco_data = Reco_calculate_columns(:,1);
Reco_data = table2array(Reco_data);

Reco_DT_data = Reco_calculate_columns(:,2);
Reco_DT_data = table2array(Reco_DT_data);

DoY_data = result_data.('DoY');
Rg_data = result_data.('Rg');

% 找到 Rg_data 中不为0的位置
non_zero_indices = Rg_data ~= 0;

% 提取对应位置的 GPP_data 和 DoY_data
GPP_DT_filtered = GPP_DT_data(non_zero_indices);
GPP_f_filtered = GPP_f_data(non_zero_indices);

Reco_filtered = Reco_data(non_zero_indices);
Reco_DT_filtered = Reco_DT_data(non_zero_indices);

NEE_filtered = NEE_data(non_zero_indices);
Tair_filtered = Tair_data(non_zero_indices);
VPD_filtered = VPD_data(non_zero_indices);

LE_filtered = LE_data(non_zero_indices);

DoY_filtered = DoY_data(non_zero_indices);

% 获取唯一的日期（1到365）
unique_days = unique(DoY_filtered);

% 初始化存储结果的数组
GPP_daily_DT_total = zeros(length(unique_days), 1);
GPP_daily_f_total = zeros(length(unique_days), 1);

Reco_daily_total = zeros(length(unique_days), 1);
Reco_daily_DT_total = zeros(length(unique_days), 1);

NEE_daily_total = zeros(length(unique_days), 1);
Tair_daily_total = zeros(length(unique_days), 1);
VPD_daily_total = zeros(length(unique_days), 1);

ET_daily_total = zeros(length(unique_days), 1);

% 循环计算每天的数据（例如，计算每天的GPP平均值）
for i = 1:length(unique_days)
    % 找到当前日期对应的索引
    current_day_indices = DoY_filtered == unique_days(i);
    
    % 提取对应的 GPP 数据
    GPP_DT_day = GPP_DT_filtered(current_day_indices);
    GPP_f_day = GPP_f_filtered(current_day_indices);

    Reco_day = Reco_filtered(current_day_indices);
    Reco_DT_day = Reco_DT_filtered(current_day_indices);

    NEE_day = NEE_filtered(current_day_indices);
    Tair_day = Tair_filtered(current_day_indices);
    VPD_day = VPD_filtered(current_day_indices);

    LE_day = LE_filtered(current_day_indices);
    LE_day(LE_day < 0) = 0;

    % 检查 Tair_day 中的负值或零值，并替换为 NaN
    Tair_day(Tair_day <= 0) = NaN;
    VPD_day(VPD_day <= 0) = NaN;
    
    % 使用 nanmean 计算 Tair_day 和 VPD_day 的均值，忽略 NaN 值
    Tair_day_mean = nanmean(Tair_day);
    VPD_day_mean = nanmean(VPD_day);

    seconds_per_halfhour = 1800;  % 每个半小时有1800秒

    % 计算每一天的GPP总量
    GPP_DT_daily = sum(GPP_DT_day) * seconds_per_halfhour;
    GPP_f_daily = sum(GPP_f_day) * seconds_per_halfhour;
    Reco_daily = sum(Reco_day) * seconds_per_halfhour;
    Reco_DT_daily = sum(Reco_DT_day) * seconds_per_halfhour;

    NEE_daily = sum(NEE_day) * seconds_per_halfhour;

    LE_daily = sum(LE_day) * seconds_per_halfhour;

    % 转换为 gC·m⁻²·day⁻¹，使用12 µg C对应每1 µmol CO₂
    GPP_DT_daily_gC = GPP_DT_daily * 12 * 10^(-6);
    GPP_f_daily_gC = GPP_f_daily * 12 * 10^(-6);
    Reco_daily_gC = Reco_daily * 12 * 10^(-6);
    Reco_DT_daily_gC = Reco_DT_daily * 12 * 10^(-6);

    NEE_daily_gC = NEE_daily * 12 * 10^(-6);

    % 参数
    lambda = 2.45e6; % 水的蒸发潜热 (J/kg)
    
    % 计算蒸散发
    ET_daily = LE_daily / lambda ; % 每半小时的蒸散发量 (mm)
    
    GPP_daily_DT_total(i) = GPP_DT_daily_gC;
    GPP_daily_f_total(i) = GPP_f_daily_gC;
    Reco_daily_total(i) = Reco_daily_gC;
    Reco_daily_DT_total(i) = Reco_DT_daily_gC;

    NEE_daily_total(i) = NEE_daily_gC;
    Tair_daily_total(i) = Tair_day_mean;
    VPD_daily_total(i) = VPD_day_mean;

    ET_daily_total(i) = ET_daily;
end

% 拼接列向量
daily_totals_matrix = [GPP_daily_f_total, GPP_daily_DT_total, Reco_daily_total, Reco_daily_DT_total, ...
                       NEE_daily_total, Tair_daily_total, VPD_daily_total, ET_daily_total, uWUE_T_data, P_data];

% 矩阵的行数
num_days = size(daily_totals_matrix, 1);

% 每8天分为一组，计算组数
num_groups = ceil(num_days / 8);

% 初始化结果矩阵
combined_result_matrix = zeros(num_groups, size(daily_totals_matrix, 2));

% 循环计算每组的均值
for group_idx = 1:num_groups
    % 确定当前组的行索引范围
    start_idx = (group_idx - 1) * 8 + 1;
    end_idx = min(group_idx * 8, num_days); % 确保索引不超出范围

    % 提取当前组的数据
    group_data = daily_totals_matrix(start_idx:end_idx, :);

    % 计算均值（忽略 NaN）
    combined_result_matrix(group_idx, :) = mean(group_data, 'omitnan');
end

% 假设 daily_totals_matrix 已经有对应的列名
column_names = {'GPP_f', 'GPP_DT', 'Reco', 'Reco_DT', 'NEE', 'Tair', 'VPD', 'ET', 'uWUE_T', 'P'};

% 转换结果矩阵为 table 并添加列名
combined_result_table = array2table(combined_result_matrix, ...
    'VariableNames', column_names);

output_path = 'N:\GZW\站点数据_1217\千烟洲';
% % 保存为 CSV 文件
% csv_filename = fullfile(output_path, sprintf('combined_result_table_%d_%s.csv', year, site_name));
% writetable(combined_result_table, csv_filename);

% 保存为 MAT 文件
save(fullfile(output_path, sprintf('combined_result_table_%d_%s_new.mat',year,site_name)), 'combined_result_table');