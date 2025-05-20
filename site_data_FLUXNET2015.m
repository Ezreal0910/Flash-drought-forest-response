% clc;
% clear;
% 
% filename = "N:\GZW\FLUXNET2015\FLX_DE-Tha_FLUXNET2015_FULLSET_HH_1996-2014_1-4.csv";
% 
% data = readtable(filename);
% 
% % 提取第35089到第52608行的数据
% extracted_data = data(35089:52608, :);
% 
% timestamp = extracted_data{:, 2};
% 
% % 初始化存储转换结果的变量
% datetime_values = NaT(length(timestamp), 1);  % 创建空的 NaT 数组
% 
% % 遍历每个时间戳并进行处理
% for i = 1:length(timestamp)
%     % 将时间戳转换为字符串，并保留三位小数
%     timestamp_str = num2str(timestamp(i), '%.7f');  % 保留7位小数（以确保处理毫秒）
%     
%     % 确保是字符类型，且去除小数点后的部分
%     timestamp_str = char(timestamp_str);
%     
%     % 将 timestamp_str 变成一个合适的日期格式
%     date_part = timestamp_str(1:8);  % 'yyyyMMdd'
%     time_part = timestamp_str(9:12); % 'HHmm'
%     millis_part = timestamp_str(13:end); % 'SSS'
%     
%     % 拼接新的时间戳格式
%     formatted_timestamp = strcat(date_part, time_part, millis_part);
%     
%     % 使用 datetime 函数将字符串转换为 datetime 格式
%     datetime_values(i) = datetime(formatted_timestamp, 'InputFormat', 'yyyyMMddHHmm.SSS');
% end
% 
% % 计算年份的第几天
% day_of_year = day(datetime_values, 'dayofyear');  % 获取该日期是年份中的第几天
% day_of_year(end) = 366;
% 
% % 获取日期部分 (去掉时间)
% midnight = datetime_values;
% midnight.Hour = 0;
% midnight.Minute = 0;
% midnight.Second = 0;
% 
% % 计算从当天午夜到每个时间戳的分钟数
% minutes_since_midnight = minutes(datetime_values - midnight);
% 
% % 计算每个时间戳是当天的第几个小时
% hour_index = minutes_since_midnight / 60;
% hour_decimal = round(hour_index, 1);

NEE_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'NEE_'));
NEE_data = NEE_columns(:,1);
NEE_data = table2array(NEE_data);
% 
% LE_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'LE_'));
% LE_data = LE_columns(:,1);
% LE_data = table2array(LE_data);
% 
% H_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'H_'));
% H_data = H_columns(:,1);
% H_data = table2array(H_data);
% 
% SW_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'SW_'));
% SW_data = SW_columns(:,2);
% SW_data = table2array(SW_data);
% 
% TS_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'TS_'));
% TS_data = TS_columns(:,1);
% TS_data = table2array(TS_data);
% 
% TA_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'TA_'));
% TA_data = TA_columns(:,1);
% TA_data = table2array(TA_data);
% 
% RH_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'RH'));
% RH_data = RH_columns(:,1);
% RH_data = table2array(RH_data);
% 
% % 初始化 VPD 数组
% VPD_calculate = zeros(size(RH_data));
% 
% % 计算 VPD
% for i = 1:length(RH_data)
%     % 获取当前温度和湿度
%     T = TA_data(i);  % 当前温度
%     rH = RH_data(i); % 当前相对湿度
% 
%     % 计算饱和蒸气压 e_s (单位：kPa)
%     e_s = 0.6108 * exp((17.27 * T) / (T + 237.3));
%     
%     % 计算实际蒸气压 e_a (单位：kPa)
%     e_a = (rH / 100) * e_s;
%     
%     % 计算 VPD (单位：kPa)
%     VPD_calculate(i) = e_s - e_a;
% end
% 
% VPD_calculate = VPD_calculate * 10;
% 
% VPD_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'VPD'));
% VPD_data = VPD_columns(:,1);
% VPD_data = table2array(VPD_data);
% 
% WS_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'WS'));
% canopy_wind_speed = WS_columns(:,1);
% canopy_wind_speed = table2array(canopy_wind_speed);
% 
% % 风速测量的高度
% z = 40;  % 风速测量高度 39.6 m
% d = 20;
% 
% % 假设粗糙度长度为 1.5 m（可以根据实际情况进行调整）
% z0 = 1.5;
% 
% % 科尔门常数
% kappa = 0.41;
% 
% % 计算摩擦速度
% USTAR_calculate = canopy_wind_speed * kappa / log(z / z0);
% USTAR_calculate(USTAR_calculate < 0) = -9999;
% 
% USTAR_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'USTAR'));
% USTAR_data = USTAR_columns(:,1);
% USTAR_data = table2array(USTAR_data);
% 
% % t = (canopy_wind_speed .* kappa) ./ USTAR_data;
% 
% 
% % 1. 获取年份
% Year = repmat(1998, length(datetime_values), 1);
% 
% % 2. 获取 DoY（年中的第几天）
% DoY = day_of_year;
% 
% % 3. 获取小时数（已计算过 hour_decimal）
% Hour = hour_decimal;
% 
% % 4. 获取其他数据列
% NEE = NEE_data;
% LE = LE_data;
% H = H_data;
% Rg = SW_data;   % 假设 SW_data 代表 Rg 数据（如果不对，调整此部分）
% Tair = TA_data;
% Tsoil = TS_data;
% rH = RH_data;
% VPD = VPD_data;  % VPD 已经计算过
% Ustar = USTAR_calculate;
% 
% % 5. 将所有数据按列合并到一个矩阵中
% output_matrix = [Year, DoY, Hour, NEE, LE, H, Rg, Tair, Tsoil, rH, VPD, Ustar];
% 
% % 单位行 (单位数组改为字符数组，以便与矩阵合并)
% units_row = {'-', '-', '-', 'umolm-2s-1', 'Wm-2', 'WWm-2', 'Wm-2', 'degC', 'degC', '%', 'hPa', 'ms-1'};
% 
% % 列名 (列名也转换为字符数组)
% headers = {'Year', 'DoY', 'Hour', 'NEE', 'LE', 'H', 'Rg', 'Tair', 'Tsoil', 'rH', 'VPD', 'Ustar'};
% 
% % 将数值矩阵转换为单元格数组
% output_cell = num2cell(output_matrix);
% 
% % 合并列名、单位行和数据
% combined_matrix = [headers; units_row; output_cell];
% 
% % 指定保存路径
% save_path = 'N:\GZW\FLUXNET2015\对照数据\output_data_with_units_Ustar2.txt';
% 
% % 保存 combined_matrix 为文本文件，使用制表符作为分隔符
% writecell(combined_matrix, save_path, 'Delimiter', '\t');

GPP_columns = extracted_data(:, contains(extracted_data.Properties.VariableNames, 'GPP'));
GPP_NT_data = GPP_columns(:,12);
GPP_NT_data = table2array(GPP_NT_data);
GPP_DT_data = GPP_columns(:,34);
GPP_DT_data = table2array(GPP_DT_data);

result_filename = "N:\GZW\FLUXNET2015\对照数据\Result3\output.txt";

% 使用 readtable 读取数据，自动处理分隔符
result_data = readtable(result_filename);

GPP_calculate_columns = result_data(:, contains(result_data.Properties.VariableNames, 'GPP'));


