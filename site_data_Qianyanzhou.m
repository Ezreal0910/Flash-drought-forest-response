clc;
clear;

% 定义文件夹路径
meteoro_data_path = 'N:\GZW\站点数据_1217\千烟洲\气象数据\30分钟';
flux_data_path = 'N:\GZW\站点数据_1217\千烟洲\通量数据\30分钟';

flux_data_path2 = "N:\GZW\站点数据_1217\千烟洲\FLX_CN-Qia_FLUXNET2015_FULLSET_2003-2005_1-4\FLX_CN-Qia_FLUXNET2015_FULLSET_HH_2003-2005_1-4.csv";
flux_data2 = readtable(flux_data_path2, 'VariableNamingRule', 'preserve');

% 提取第一列数据
first_column = flux_data2{:, 1}; % 如果是字符串或数值列

% 将第一列数据转换为字符串，并提取前四位
if isnumeric(first_column)
    first_column = string(first_column); % 如果是数值列，转换为字符串
end
first_four_chars = extractBefore(first_column, 5); % 提取前四位字符

% 使用 findgroups 创建分组索引
group_indices = findgroups(first_four_chars);
column_names = flux_data2.Properties.VariableNames;

grouped_data = splitapply(@(varargin) {array2table([varargin{:}], 'VariableNames', column_names)}, ...
                          table2array(flux_data2), group_indices);

output_path = 'N:\GZW\站点数据_1217\千烟洲';

% 获取气象数据文件夹中的所有文件
meteoro_files = dir(fullfile(meteoro_data_path, '*'));

% 获取通量数据文件夹中的所有文件
flux_files = dir(fullfile(flux_data_path, '*'));

% 过滤掉目录（只保留文件）
meteoro_files = meteoro_files(~[meteoro_files.isdir]);  % 只保留文件，不包括文件夹
flux_files = flux_files(~[flux_files.isdir]);  % 只保留文件，不包括文件夹

for i = 2:length(meteoro_files)
    % 获取文件的完整路径
    meteoro_file_path = fullfile(meteoro_data_path, meteoro_files(i).name);
    flux_file_path = fullfile(flux_data_path, flux_files(i).name);

    target_data = grouped_data{i};

    Rn_data = target_data.('NETRAD');
    H_data = target_data.('H_F_MDS');
    
    LE_data = Rn_data - H_data;
   
    % 读取 CSV 文件，并保留原始列标题
    meteoro_data = readtable(meteoro_file_path, 'VariableNamingRule', 'preserve');
    flux_data = readtable(flux_file_path, 'VariableNamingRule', 'preserve');

    year_col = meteoro_data.('年');
    month_col = meteoro_data.('月');
    day_col = meteoro_data.('日');
    hour_col = meteoro_data.('时');
    minute_col = meteoro_data.('分');
    second_col = meteoro_data.('秒');

    % 转换为 datetime 数组
    datetime_array = datetime(year_col, month_col, day_col, ...
        hour_col, minute_col, second_col);

    % 计算年份的第几天
    day_of_year = day(datetime_array, 'dayofyear');  % 获取该日期是年份中的第几天
    day_of_year(end) = day_of_year(end-1);
    
    % 获取日期部分 (去掉时间)
    midnight = datetime_array;
    midnight.Hour = 0;
    midnight.Minute = 0;
    midnight.Second = 0;
    
    % 计算从当天午夜到每个时间戳的分钟数
    minutes_since_midnight = minutes(datetime_array - midnight);
    
    % 计算每个时间戳是当天的第几个小时
    hour_index = minutes_since_midnight / 60;
    hour_decimal = round(hour_index, 1);
    hour_decimal(end) = 24;

    NEE_data = flux_data.('NEE');
    NEE_data = NEE_data * 22.7;
%     LE_data = flux_data.('LE');
    Hs_data = flux_data.('Hs');

    TS_data = meteoro_data.('一层土壤温度');
    TA_data = meteoro_data.('冠层上方空气温度');
    RH_data = meteoro_data.('冠层上方空气湿度');
    SW_data = meteoro_data.('太阳辐射');

    % 初始化 VPD 数组
    VPD_calculate = zeros(size(RH_data));
    
    % 计算 VPD
    for j = 1:length(RH_data)
        % 获取当前温度和湿度
        T = TA_data(j);  % 当前温度
        rH = RH_data(j); % 当前相对湿度
    
        % 计算饱和蒸气压 e_s (单位：kPa)
        e_s = 0.6108 * exp((17.27 * T) / (T + 237.3));
        
        % 计算实际蒸气压 e_a (单位：kPa)
        e_a = (rH / 100) * e_s;
        
        % 计算 VPD (单位：kPa)
        VPD_calculate(j) = e_s - e_a;
    end
    
    VPD_calculate = VPD_calculate * 10;

    WS_data = meteoro_data.('冠层上方风速');

    % 风速测量的高度
    z = 39.6;  % 风速测量高度 39.6 m
    
    % 假设粗糙度长度为 1.5 m（可以根据实际情况进行调整）
    z0 = 1.5;
    
    % 科尔门常数
    kappa = 0.41;
    
    % 计算摩擦速度
    USTAR_calculate = WS_data * kappa / log(z / z0);
    USTAR_calculate(USTAR_calculate < 0) = -9999;

    % 1. 获取年份
    Year = meteoro_data.("年");
    
    % 2. 获取 DoY（年中的第几天）
    DoY = day_of_year;
    
    % 3. 获取小时数（已计算过 hour_decimal）
    Hour = hour_decimal;
    
    % 4. 获取其他数据列
    NEE = NEE_data;
    LE = LE_data;
    H = Hs_data;
    Rg = SW_data;   % 假设 SW_data 代表 Rg 数据（如果不对，调整此部分）
    Tair = TA_data;
    Tsoil = TS_data;
    rH = RH_data;
    VPD = VPD_calculate;  % VPD 已经计算过
    Ustar = USTAR_calculate;

    % 对每个变量四舍五入到4位小数
    Year = round(Year, 4);
    DoY = round(DoY, 4);
    Hour = round(Hour, 4);
    NEE = round(NEE, 4);
    LE = round(LE, 4);
    H = round(H, 4);
    Rg = round(Rg, 4);
    Tair = round(Tair, 4);
    Tsoil = round(Tsoil, 4);
    rH = round(rH, 4);
    VPD = round(VPD, 4);
    Ustar = round(Ustar, 4);
    
    % 5. 将所有数据按列合并到一个矩阵中
    output_matrix = [Year, DoY, Hour, NEE, LE, H, Rg, Tair, Tsoil, rH, VPD, Ustar];
    
    % 单位行 (单位数组改为字符数组，以便与矩阵合并)
    units_row = {'-', '-', '-', 'umolm-2s-1', 'Wm-2', 'Wm-2', 'Wm-2', 'degC', 'degC', '%', 'hPa', 'ms-1'};
    
    % 列名 (列名也转换为字符数组)
    headers = {'Year', 'DoY', 'Hour', 'NEE', 'LE', 'H', 'Rg', 'Tair', 'Tsoil', 'rH', 'VPD', 'Ustar'};
    
    % 将数值矩阵转换为单元格数组
    output_cell = num2cell(output_matrix);
    
    % 合并列名、单位行和数据
    combined_matrix = [headers; units_row; output_cell];

    % 构建新保存路径，假设你想将文件保存在 `Result` 文件夹下
    save_path = fullfile(output_path, 'Result');  % 在原文件夹下创建一个新文件夹
    
    % 如果文件夹不存在，则创建它
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end

    full_path = fullfile(save_path, sprintf('%d_output_data_new.txt', Year(1)));
    
    % 保存 combined_matrix 为文本文件，使用制表符作为分隔符
    writecell(combined_matrix, full_path, 'Delimiter', '\t');
end
