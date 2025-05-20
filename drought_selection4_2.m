clc;
clear;

%   分类为季节的同时，也统计了干旱事件的开始和结束时间
%   还可以统计  人工林  和   自然林

% 定义季节区间
interval_spring = [3, 4, 5];
interval_summer = [6, 7, 8];
interval_autumn = [9, 10, 11];
interval_winter = [12, 1, 2];

load("G:\GZW\Drought_response\date\dates.mat");

% 初始化一个 920x3 的二维矩阵，用于存储年、月和日
dateMatrix = zeros(920, 3);

% 使用循环将字符数组中的每个日期字符串转换为二维矩阵的行
for i = 1:920
    % 从字符数组中提取年、月和日部分
    yearStr = dates(i, 1:4);
    monthStr = dates(i, 6:7);
    dayStr = dates(i, 9:10);
    
    % 将提取的字符串转换为数字，并存储在二维矩阵中
    dateMatrix(i, 1) = str2double(yearStr);
    dateMatrix(i, 2) = str2double(monthStr);
    dateMatrix(i, 3) = str2double(dayStr);
end

tif_FileName = "G:\GZW\Landcover\Landcover_forest.tif";
forest_data = imread(tif_FileName);

% tif_FileName = "E:\GZW\Natural.tif";
% forest_data2 = imread(tif_FileName);
% 
% tif_FileName = "E:\GZW\Planted.tif";
% forest_data3 = imread(tif_FileName);

forest_data(forest_data == -128) = 0;

% 获取 FD_events
FD_events_path = "F:\GZW\Output\Drought_events\FD_events";
FD_events_files = dir(fullfile(FD_events_path, '*.mat'));

% 获取 SD_events
SD_events_path = "F:\GZW\Output\Drought_events\SD_events";
SD_events_files = dir(fullfile(SD_events_path, '*.mat'));

% 获取 ESI 文件
ESI_folder_path = "F:\GZW\Output\ESI_Assign";
ESI_files = dir(fullfile(ESI_folder_path, '*.mat'));

% 获取 ΔESI 文件
DESI_folder_path = "F:\GZW\Output\ΔESI";
DESI_files = dir(fullfile(DESI_folder_path, '*.mat'));

cols_per_block = 620;

result_FD_Events_spring = [];
result_FD_Events_spring_start = []; 
result_FD_Events_spring_end = []; 

result_FD_Events_summer = [];
result_FD_Events_summer_start = [];
result_FD_Events_summer_end = [];

result_FD_Events_autumn = [];
result_FD_Events_autumn_start = [];
result_FD_Events_autumn_end = [];

result_FD_Events_winter = [];
result_FD_Events_winter_start = [];
result_FD_Events_winter_end = [];

result_SD_Events_spring = [];
result_SD_Events_spring_start = [];
result_SD_Events_spring_end = [];

result_SD_Events_summer = [];
result_SD_Events_summer_start = [];
result_SD_Events_summer_end = [];

result_SD_Events_autumn = [];
result_SD_Events_autumn_start = [];
result_SD_Events_autumn_end = [];

result_SD_Events_winter = [];
result_SD_Events_winter_start = [];
result_SD_Events_winter_end = [];

result_Flash_Drought_Events_spring = {};
result_Flash_Drought_Events_spring_start = {};
result_Flash_Drought_Events_spring_end = {};

result_Flash_Drought_Events_summer = {};
result_Flash_Drought_Events_summer_start = {};
result_Flash_Drought_Events_summer_end = {};

result_Flash_Drought_Events_autumn = {};
result_Flash_Drought_Events_autumn_start = {};
result_Flash_Drought_Events_autumn_end = {};

result_Flash_Drought_Events_winter = {};
result_Flash_Drought_Events_winter_start = {};
result_Flash_Drought_Events_winter_end = {};

result_Slow_Drought_Events_spring = {};
result_Slow_Drought_Events_spring_start = {};
result_Slow_Drought_Events_spring_end = {};

result_Slow_Drought_Events_summer = {};
result_Slow_Drought_Events_summer_start = {};
result_Slow_Drought_Events_summer_end = {};

result_Slow_Drought_Events_autumn = {};
result_Slow_Drought_Events_autumn_start = {};
result_Slow_Drought_Events_autumn_end = {};

result_Slow_Drought_Events_winter = {};
result_Slow_Drought_Events_winter_start = {};
result_Slow_Drought_Events_winter_end = {};

result_FD_Events_slow_spring = [];
result_FD_Events_slow_spring_start = [];
result_FD_Events_slow_spring_end = [];

result_FD_Events_slow_summer = [];
result_FD_Events_slow_summer_start = [];
result_FD_Events_slow_summer_end = [];

result_FD_Events_slow_autumn = [];
result_FD_Events_slow_autumn_start = [];
result_FD_Events_slow_autumn_end = [];

result_FD_Events_slow_winter = [];
result_FD_Events_slow_winter_start = [];
result_FD_Events_slow_winter_end = [];

result_SD_Events_slow_spring = [];
result_SD_Events_slow_spring_start = [];
result_SD_Events_slow_spring_end = [];

result_SD_Events_slow_summer = [];
result_SD_Events_slow_summer_start = [];
result_SD_Events_slow_summer_end = [];

result_SD_Events_slow_autumn = [];
result_SD_Events_slow_autumn_start = [];
result_SD_Events_slow_autumn_end = [];

result_SD_Events_slow_winter = [];
result_SD_Events_slow_winter_start = [];
result_SD_Events_slow_winter_end = [];

% 遍历每个文件
for t = 8 : 10

    start_col = (t-1)*cols_per_block + 1;
    end_col = min(start_col+cols_per_block-1 , 6200);
    
    col = end_col - start_col + 1;

    load(fullfile(FD_events_path, FD_events_files(t).name));
    load(fullfile(SD_events_path, SD_events_files(t).name));

    forest_extract = forest_data(:, start_col:end_col);
    
    % 使用 all 函数检查整个矩阵是否都为零
    if all(forest_extract(:) == 0)
        % 如果整个矩阵都为空，跳过本次循环
        disp('整个矩阵为空，跳过循环');
        [row_count, col_count] = size(forest_extract);
        leaped_Events1 = zeros(row_count, col_count);
        leaped_Events2 = cell(row_count, col_count);

        result_FD_Events_spring = [result_FD_Events_spring, leaped_Events1];
        result_FD_Events_spring_start = [result_FD_Events_spring_start, leaped_Events1];
        result_FD_Events_spring_end = [result_FD_Events_spring_end, leaped_Events1];

        result_SD_Events_spring = [result_SD_Events_spring, leaped_Events1];
        result_SD_Events_spring_start = [result_SD_Events_spring_start, leaped_Events1];
        result_SD_Events_spring_end = [result_SD_Events_spring_end, leaped_Events1];

        result_FD_Events_summer = [result_FD_Events_summer, leaped_Events1];
        result_FD_Events_summer_start = [result_FD_Events_summer_start, leaped_Events1];
        result_FD_Events_summer_end = [result_FD_Events_summer_end, leaped_Events1];

        result_SD_Events_summer = [result_SD_Events_summer, leaped_Events1];
        result_SD_Events_summer_start = [result_SD_Events_summer_start, leaped_Events1];
        result_SD_Events_summer_end = [result_SD_Events_summer_end, leaped_Events1];

        result_FD_Events_autumn = [result_FD_Events_autumn, leaped_Events1];
        result_FD_Events_autumn_start = [result_FD_Events_autumn_start, leaped_Events1];
        result_FD_Events_autumn_end = [result_FD_Events_autumn_end, leaped_Events1];

        result_SD_Events_autumn = [result_SD_Events_autumn, leaped_Events1];
        result_SD_Events_autumn_start = [result_SD_Events_autumn_start, leaped_Events1];
        result_SD_Events_autumn_end = [result_SD_Events_autumn_end, leaped_Events1];

        result_FD_Events_winter = [result_FD_Events_winter, leaped_Events1];
        result_FD_Events_winter_start = [result_FD_Events_winter_start, leaped_Events1];
        result_FD_Events_winter_end = [result_FD_Events_winter_end, leaped_Events1];

        result_SD_Events_winter = [result_SD_Events_winter, leaped_Events1];
        result_SD_Events_winter_start = [result_SD_Events_winter_start, leaped_Events1];
        result_SD_Events_winter_end = [result_SD_Events_winter_end, leaped_Events1];

        result_Flash_Drought_Events_spring = [result_Flash_Drought_Events_spring, leaped_Events2];
        result_Flash_Drought_Events_spring_start = [result_Flash_Drought_Events_spring_start, leaped_Events2];
        result_Flash_Drought_Events_spring_end = [result_Flash_Drought_Events_spring_end, leaped_Events2];

        result_Slow_Drought_Events_spring = [result_Slow_Drought_Events_spring, leaped_Events2];
        result_Slow_Drought_Events_spring_start = [result_Slow_Drought_Events_spring_start, leaped_Events2];
        result_Slow_Drought_Events_spring_end = [result_Slow_Drought_Events_spring_end, leaped_Events2];

        result_Flash_Drought_Events_summer = [result_Flash_Drought_Events_summer, leaped_Events2];
        result_Flash_Drought_Events_summer_start = [result_Flash_Drought_Events_summer_start, leaped_Events2];
        result_Flash_Drought_Events_summer_end = [result_Flash_Drought_Events_summer_end, leaped_Events2];

        result_Slow_Drought_Events_summer = [result_Slow_Drought_Events_summer, leaped_Events2];
        result_Slow_Drought_Events_summer_start = [result_Slow_Drought_Events_summer_start, leaped_Events2];
        result_Slow_Drought_Events_summer_end = [result_Slow_Drought_Events_summer_end, leaped_Events2];

        result_Flash_Drought_Events_autumn = [result_Flash_Drought_Events_autumn, leaped_Events2];
        result_Flash_Drought_Events_autumn_start = [result_Flash_Drought_Events_autumn_start, leaped_Events2];
        result_Flash_Drought_Events_autumn_end = [result_Flash_Drought_Events_autumn_end, leaped_Events2];

        result_Slow_Drought_Events_autumn = [result_Slow_Drought_Events_autumn, leaped_Events2];
        result_Slow_Drought_Events_autumn_start = [result_Slow_Drought_Events_autumn_start, leaped_Events2];
        result_Slow_Drought_Events_autumn_end = [result_Slow_Drought_Events_autumn_end, leaped_Events2];

        result_Flash_Drought_Events_winter = [result_Flash_Drought_Events_winter, leaped_Events2];
        result_Flash_Drought_Events_winter_start = [result_Flash_Drought_Events_winter_start, leaped_Events2];
        result_Flash_Drought_Events_winter_end = [result_Flash_Drought_Events_winter_end, leaped_Events2];

        result_Slow_Drought_Events_winter = [result_Slow_Drought_Events_winter, leaped_Events2];
        result_Slow_Drought_Events_winter_start = [result_Slow_Drought_Events_winter_start, leaped_Events2];
        result_Slow_Drought_Events_winter_end = [result_Slow_Drought_Events_winter_end, leaped_Events2];

        result_FD_Events_slow_spring = [result_FD_Events_slow_spring, leaped_Events1];
        result_FD_Events_slow_spring_start = [result_FD_Events_slow_spring_start, leaped_Events1];
        result_FD_Events_slow_spring_end = [result_FD_Events_slow_spring_end, leaped_Events1];

        result_SD_Events_slow_spring = [result_SD_Events_slow_spring, leaped_Events1];
        result_SD_Events_slow_spring_start = [result_SD_Events_slow_spring_start, leaped_Events1];
        result_SD_Events_slow_spring_end = [result_SD_Events_slow_spring_end, leaped_Events1];

        result_FD_Events_slow_summer = [result_FD_Events_slow_summer, leaped_Events1];
        result_FD_Events_slow_summer_start = [result_FD_Events_slow_summer_start, leaped_Events1];
        result_FD_Events_slow_summer_end = [result_FD_Events_slow_summer_end, leaped_Events1];

        result_SD_Events_slow_summer = [result_SD_Events_slow_summer, leaped_Events1];
        result_SD_Events_slow_summer_start = [result_SD_Events_slow_summer_start, leaped_Events1];
        result_SD_Events_slow_summer_end = [result_SD_Events_slow_summer_end, leaped_Events1];

        result_FD_Events_slow_autumn = [result_FD_Events_slow_autumn, leaped_Events1];
        result_FD_Events_slow_autumn_start = [result_FD_Events_slow_autumn_start, leaped_Events1];
        result_FD_Events_slow_autumn_end = [result_FD_Events_slow_autumn_end, leaped_Events1];

        result_SD_Events_slow_autumn = [result_SD_Events_slow_autumn, leaped_Events1];
        result_SD_Events_slow_autumn_start = [result_SD_Events_slow_autumn_start, leaped_Events1];
        result_SD_Events_slow_autumn_end = [result_SD_Events_slow_autumn_end, leaped_Events1];

        result_FD_Events_slow_winter = [result_FD_Events_slow_winter, leaped_Events1];
        result_FD_Events_slow_winter_start = [result_FD_Events_slow_winter_start, leaped_Events1];
        result_FD_Events_slow_winter_end = [result_FD_Events_slow_winter_end, leaped_Events1];

        result_SD_Events_slow_winter = [result_SD_Events_slow_winter, leaped_Events1];
        result_SD_Events_slow_winter_start = [result_SD_Events_slow_winter_start, leaped_Events1];
        result_SD_Events_slow_winter_end = [result_SD_Events_slow_winter_end, leaped_Events1];
        continue
    end

    all_DESI = [];

    all_ESI = [];

    for i = 1 : numel(ESI_files)

        load(fullfile(DESI_folder_path, DESI_files(i).name));
        
        load(fullfile(ESI_folder_path, ESI_files(i).name));

        % 获取对应列的ΔESI
        DESI_3D = DESI_data(:, start_col:end_col, :);

        all_DESI = cat(3, all_DESI, DESI_3D);

        % 获取对应列的ESI
        ESI_3D = ESI_data(:, start_col:end_col, :);

        all_ESI = cat(3, all_ESI, ESI_3D);
    end
    
    DESI_2D = reshape(all_DESI, [], size(all_DESI, 3));
    ESI_2D = reshape(all_ESI, [], size(all_ESI, 3));

    % 找到不为 0 的位置
    positions = find(forest_extract ~= 0);

    % 使用之前获取的索引从 FD_components 中提取对应位置的值
    extracted_values_FD = FD_components(positions);

    % 使用之前获取的索引从 SD_components 中提取对应位置的值
    extracted_values_SD = SD_components(positions);
    
    % 使用相同的索引从 ΔESI 中提取相应位置的所有列
    extracted_rows_DESI = DESI_2D(positions, :);

    % 使用相同的索引从 ESI 中提取相应位置的所有列
    extracted_rows_ESI = ESI_2D(positions, :);

    selected_FD_event_point_spring = zeros(3900*col,1);
    selected_FD_event_point_spring_start = zeros(3900*col,1);
    selected_FD_event_point_spring_end = zeros(3900*col,1);

    selected_FD_event_point_summer = zeros(3900*col,1);
    selected_FD_event_point_summer_start = zeros(3900*col,1);
    selected_FD_event_point_summer_end = zeros(3900*col,1);

    selected_FD_event_point_autumn = zeros(3900*col,1);
    selected_FD_event_point_autumn_start = zeros(3900*col,1);
    selected_FD_event_point_autumn_end = zeros(3900*col,1);

    selected_FD_event_point_winter = zeros(3900*col,1);
    selected_FD_event_point_winter_start = zeros(3900*col,1);
    selected_FD_event_point_winter_end = zeros(3900*col,1);

    selected_FD_event_point_slow_spring = zeros(3900*col,1);
    selected_FD_event_point_slow_spring_start = zeros(3900*col,1);
    selected_FD_event_point_slow_spring_end = zeros(3900*col,1);

    selected_FD_event_point_slow_summer = zeros(3900*col,1);
    selected_FD_event_point_slow_summer_start = zeros(3900*col,1);
    selected_FD_event_point_slow_summer_end = zeros(3900*col,1);

    selected_FD_event_point_slow_autumn = zeros(3900*col,1);
    selected_FD_event_point_slow_autumn_start = zeros(3900*col,1);
    selected_FD_event_point_slow_autumn_end = zeros(3900*col,1);

    selected_FD_event_point_slow_winter = zeros(3900*col,1);
    selected_FD_event_point_slow_winter_start = zeros(3900*col,1);
    selected_FD_event_point_slow_winter_end = zeros(3900*col,1);

    selected_many_FD_event_spring = cell(3900*col,1);
    selected_many_FD_event_spring_start = cell(3900*col,1);
    selected_many_FD_event_spring_end = cell(3900*col,1);

    selected_many_FD_event_summer = cell(3900*col,1);
    selected_many_FD_event_summer_start = cell(3900*col,1);
    selected_many_FD_event_summer_end = cell(3900*col,1);

    selected_many_FD_event_autumn = cell(3900*col,1);
    selected_many_FD_event_autumn_start = cell(3900*col,1);
    selected_many_FD_event_autumn_end = cell(3900*col,1);

    selected_many_FD_event_winter = cell(3900*col,1);
    selected_many_FD_event_winter_start = cell(3900*col,1);
    selected_many_FD_event_winter_end = cell(3900*col,1);
    
    for m = 1:numel(extracted_values_FD)

        FD_Events = extracted_values_FD{m};
        SD_Events = extracted_values_SD{m};

        if isempty(FD_Events)
            continue;
        end

        % 获取 NumObjects 字段的值
        FD_numObjects = FD_Events.NumObjects;

        if FD_numObjects == 0
            % 如果FD_numObjects为0，跳出循环
            continue;
        end
        
        DESI_contrast_values = extracted_rows_DESI(m, :);
        ESI_contrast_values = extracted_rows_ESI(m, :);

        position_index = positions(m);

        selected_FD_Event = [];  % 创建一个空的矩阵
        
        for e1 = 1:numel(FD_Events.PixelIdxList)
            FD_Event = FD_Events.PixelIdxList{e1};

            if isempty(FD_Event)
                continue
            end

            FD_Event_ESI_values = ESI_contrast_values(FD_Event);
            [FD_min_value, FD_min_index] = min(FD_Event_ESI_values);

            FD_Event_min_ESI_period = FD_Event(FD_min_index);
            FD_Event_first_ESI_period = FD_Event(1);

            if FD_Event_min_ESI_period <= 12 || FD_Event_min_ESI_period > 908
                % 最小值的周期在两端，则跳过当前循环
                continue;
            end

            FD_ESI_decrease_idx = FD_Event(1:FD_min_index);
            FD_DESI_decrease_idx = FD_ESI_decrease_idx(1:end-1);
            FD_Event_DESI_values = DESI_contrast_values(FD_DESI_decrease_idx);
            mean_FD_Event_DESI_values = mean(FD_Event_DESI_values);

            % 获取当前元素的第一个和最后一个值
            FD_Event_drought_start = FD_Event(1);
            FD_Event_drought_end = FD_Event(end);

            % 将这两个值存储为当前行
            current_row = [FD_Event_first_ESI_period, FD_Event_min_ESI_period, FD_min_value, mean_FD_Event_DESI_values, FD_Event_drought_start, FD_Event_drought_end];

            selected_FD_Event = [selected_FD_Event; current_row];
        end

        if isempty(selected_FD_Event)
            continue;
        end
        
        first_column_spring = [];
        first_column_spring_start = [];
        first_column_spring_end = [];

        first_column_summer = [];
        first_column_summer_start = [];
        first_column_summer_end = [];

        first_column_autumn = [];
        first_column_autumn_start = [];
        first_column_autumn_end = [];

        first_column_winter = [];
        first_column_winter_start = [];
        first_column_winter_end = [];

        % 记录添加的次数
        adding_count_spring = 0;
        adding_count_summer = 0;
        adding_count_autumn = 0;
        adding_count_winter = 0;

        % 按照第四列从小到大排列
        sorted_FD_Event = sortrows(selected_FD_Event, 4);

        % 循环遍历排序后的数组
        for f = 1:size(sorted_FD_Event, 1)
            % 获取第 i 行的数据
            drought_start = sorted_FD_Event(f, 1);
            drought_date = sorted_FD_Event(f, 2);

            % 获得前后 12 个索引号
            startIndex =  drought_date - 12;
            endIndex =  drought_date + 12;
            indices = startIndex:endIndex;

            overlapping_count = 0; % 记录重叠的次数

            for e2 = 1:numel(FD_Events.PixelIdxList)
                FD_Event = FD_Events.PixelIdxList{e2};
                
                % 检查FD_Event是否与indices有重叠
                if ~isempty(intersect(FD_Event, indices))
                    % 有重叠，增加重叠次数
                    overlapping_count = overlapping_count + 1;
                    
                    % 将该FD_Event添加到overlapping_FD_Events
                    only_one_overlapping = FD_Event;

                    % 如果重叠次数大于一个，直接进行下一次循环
                    if overlapping_count > 1
                        break;
                    end
                end
            end

            for e3 = 1:numel(SD_Events.PixelIdxList)

                SD_Event = SD_Events.PixelIdxList{e3};
                
                % 检查SD_Event是否与indices有重叠
                if ~isempty(intersect(SD_Event, indices))
                    % 有重叠，增加重叠次数
                    overlapping_count = overlapping_count + 1;

                    % 如果重叠次数大于一个，直接进行下一次循环
                    if overlapping_count > 1
                        break;
                    end
                end
            end

            % 如果重叠次数等于1
            if overlapping_count == 1
                % 检查drought_date是否包含在FD_Event中
                if ismember(drought_date, only_one_overlapping)
                    % 提取最小值对应的 FD_Event_min_ESI_period 值
                    result_FD_Event = drought_date;
                    
                    % 获取开始的月份和年份
                    start_year = dateMatrix(drought_start, 1);
                    start_month = dateMatrix(drought_start, 2);
        
                    % 获取最严重时的月份和年份
                    end_year = dateMatrix(result_FD_Event, 1);
                    end_month = dateMatrix(result_FD_Event, 2);
                    
                    % 获取开始的季节和最严重时的季节
                    start_season = findSeason(start_month, interval_spring, interval_summer, interval_autumn, interval_winter);
                    end_season = findSeason(end_month, interval_spring, interval_summer, interval_autumn, interval_winter);

                    % 干旱开始
                    FD_drought_start = sorted_FD_Event(f, 5);
                    % 干旱结束
                    FD_drought_end = sorted_FD_Event(f, 6);
                    
                    % 根据不同季节执行相应操作
                    if isequal(start_season, end_season)
                        % 在相同的季节区间内
                        switch start_season
                            case 'spring'
                                if adding_count_spring == 0
                                    selected_FD_event_point_spring(position_index) = result_FD_Event;
                                    selected_FD_event_point_spring_start(position_index) = FD_drought_start;
                                    selected_FD_event_point_spring_end(position_index) = FD_drought_end;
                                    adding_count_spring = adding_count_spring + 1;
                                end
                                % 记录多个符合条件的事件
                                first_column_spring(end+1) = result_FD_Event;
                                first_column_spring_start(end+1) = FD_drought_start;
                                first_column_spring_end(end+1) = FD_drought_end;
                            case 'summer'
                                if adding_count_summer == 0
                                    selected_FD_event_point_summer(position_index) = result_FD_Event;
                                    selected_FD_event_point_summer_start(position_index) = FD_drought_start;
                                    selected_FD_event_point_summer_end(position_index) = FD_drought_end;
                                    adding_count_summer = adding_count_summer + 1;
                                end
                                % 记录多个符合条件的事件
                                first_column_summer(end+1) = result_FD_Event;
                                first_column_summer_start(end+1) = FD_drought_start;
                                first_column_summer_end(end+1) = FD_drought_end;
                            case 'autumn'
                                if adding_count_autumn == 0
                                    selected_FD_event_point_autumn(position_index) = result_FD_Event;
                                    selected_FD_event_point_autumn_start(position_index) = FD_drought_start;
                                    selected_FD_event_point_autumn_end(position_index) = FD_drought_end;
                                    adding_count_autumn = adding_count_autumn + 1;
                                end
                                % 记录多个符合条件的事件
                                first_column_autumn(end+1) = result_FD_Event;
                                first_column_autumn_start(end+1) = FD_drought_start;
                                first_column_autumn_end(end+1) = FD_drought_end;
                            case 'winter'
                                if adding_count_winter == 0
                                    selected_FD_event_point_winter(position_index) = result_FD_Event;
                                    selected_FD_event_point_winter_start(position_index) = FD_drought_start;
                                    selected_FD_event_point_winter_end(position_index) = FD_drought_end;
                                    adding_count_winter = adding_count_winter + 1;
                                end
                                % 记录多个符合条件的事件
                                first_column_winter(end+1) = result_FD_Event;
                                first_column_winter_start(end+1) = FD_drought_start;
                                first_column_winter_end(end+1) = FD_drought_end;
                            otherwise
                                disp('Unknown season.');
                                % 添加其他处理逻辑
                        end
                    else
                        % 不在相同的季节区间内
                        disp('Invalid season intervals for start_month and/or end_month.');
                        % 添加其他处理逻辑
                    end
                end
            end
        end
        
        if ~isempty(first_column_spring)
            FD_slow_value_spring = first_column_spring(end);
            FD_slow_value_spring_start = first_column_spring_start(end);
            FD_slow_value_spring_end = first_column_spring_end(end);

            selected_FD_event_point_slow_spring(position_index) = FD_slow_value_spring;
            selected_FD_event_point_slow_spring_start(position_index) = FD_slow_value_spring_start;
            selected_FD_event_point_slow_spring_end(position_index) = FD_slow_value_spring_end;

            selected_many_FD_event_spring{position_index} = first_column_spring;
            selected_many_FD_event_spring_start{position_index} = first_column_spring_start;
            selected_many_FD_event_spring_end{position_index} = first_column_spring_end;
        end

        if ~isempty(first_column_summer)
            FD_slow_value_summer = first_column_summer(end);
            FD_slow_value_summer_start = first_column_summer_start(end);
            FD_slow_value_summer_end = first_column_summer_end(end);

            selected_FD_event_point_slow_summer(position_index) = FD_slow_value_summer;
            selected_FD_event_point_slow_summer_start(position_index) = FD_slow_value_summer_start;
            selected_FD_event_point_slow_summer_end(position_index) = FD_slow_value_summer_end;

            selected_many_FD_event_summer{position_index} = first_column_summer;
            selected_many_FD_event_summer_start{position_index} = first_column_summer_start;
            selected_many_FD_event_summer_end{position_index} = first_column_summer_end;
        end

        if ~isempty(first_column_autumn)
            FD_slow_value_autumn = first_column_autumn(end);
            FD_slow_value_autumn_start = first_column_autumn_start(end);
            FD_slow_value_autumn_end = first_column_autumn_end(end);

            selected_FD_event_point_slow_autumn(position_index) = FD_slow_value_autumn;
            selected_FD_event_point_slow_autumn_start(position_index) = FD_slow_value_autumn_start;
            selected_FD_event_point_slow_autumn_end(position_index) = FD_slow_value_autumn_end;

            selected_many_FD_event_autumn{position_index} = first_column_autumn;
            selected_many_FD_event_autumn_start{position_index} = first_column_autumn_start;
            selected_many_FD_event_autumn_end{position_index} = first_column_autumn_end;
        end

        if ~isempty(first_column_winter)
            FD_slow_value_winter = first_column_winter(end);
            FD_slow_value_winter_start = first_column_winter_start(end);
            FD_slow_value_winter_end = first_column_winter_end(end);

            selected_FD_event_point_slow_winter(position_index) = FD_slow_value_winter;
            selected_FD_event_point_slow_winter_start(position_index) = FD_slow_value_winter_start;
            selected_FD_event_point_slow_winter_end(position_index) = FD_slow_value_winter_end;

            selected_many_FD_event_winter{position_index} = first_column_winter;
            selected_many_FD_event_winter_start{position_index} = first_column_winter_start;
            selected_many_FD_event_winter_end{position_index} = first_column_winter_end;
        end
    end

    % 最快的
    % spring
    selected_FD_event_point_2D_spring = reshape(selected_FD_event_point_spring, 3900, col);
    selected_FD_event_point_2D_spring_start = reshape(selected_FD_event_point_spring_start, 3900, col);
    selected_FD_event_point_2D_spring_end = reshape(selected_FD_event_point_spring_end, 3900, col);

    result_FD_Events_spring = [result_FD_Events_spring, selected_FD_event_point_2D_spring];
    result_FD_Events_spring_start = [result_FD_Events_spring_start, selected_FD_event_point_2D_spring_start];
    result_FD_Events_spring_end = [result_FD_Events_spring_end, selected_FD_event_point_2D_spring_end];
    % summer
    selected_FD_event_point_2D_summer = reshape(selected_FD_event_point_summer, 3900, col);
    selected_FD_event_point_2D_summer_start = reshape(selected_FD_event_point_summer_start, 3900, col);
    selected_FD_event_point_2D_summer_end = reshape(selected_FD_event_point_summer_end, 3900, col);

    result_FD_Events_summer = [result_FD_Events_summer, selected_FD_event_point_2D_summer];
    result_FD_Events_summer_start = [result_FD_Events_summer_start, selected_FD_event_point_2D_summer_start];
    result_FD_Events_summer_end = [result_FD_Events_summer_end, selected_FD_event_point_2D_summer_end];
    % autumn
    selected_FD_event_point_2D_autumn = reshape(selected_FD_event_point_autumn, 3900, col);
    selected_FD_event_point_2D_autumn_start = reshape(selected_FD_event_point_autumn_start, 3900, col);
    selected_FD_event_point_2D_autumn_end = reshape(selected_FD_event_point_autumn_end, 3900, col);

    result_FD_Events_autumn = [result_FD_Events_autumn, selected_FD_event_point_2D_autumn];
    result_FD_Events_autumn_start = [result_FD_Events_autumn_start, selected_FD_event_point_2D_autumn_start];
    result_FD_Events_autumn_end = [result_FD_Events_autumn_end, selected_FD_event_point_2D_autumn_end];
    % winter
    selected_FD_event_point_2D_winter = reshape(selected_FD_event_point_winter, 3900, col);
    selected_FD_event_point_2D_winter_start = reshape(selected_FD_event_point_winter_start, 3900, col);
    selected_FD_event_point_2D_winter_end = reshape(selected_FD_event_point_winter_end, 3900, col);

    result_FD_Events_winter = [result_FD_Events_winter, selected_FD_event_point_2D_winter];
    result_FD_Events_winter_start = [result_FD_Events_winter_start, selected_FD_event_point_2D_winter_start];
    result_FD_Events_winter_end = [result_FD_Events_winter_end, selected_FD_event_point_2D_winter_end];

    % 多个
    % spring
    selected_many_FD_event_2D_spring = reshape(selected_many_FD_event_spring, 3900, col);
    selected_many_FD_event_2D_spring_start = reshape(selected_many_FD_event_spring_start, 3900, col);
    selected_many_FD_event_2D_spring_end = reshape(selected_many_FD_event_spring_end, 3900, col);

    result_Flash_Drought_Events_spring = [result_Flash_Drought_Events_spring, selected_many_FD_event_2D_spring];
    result_Flash_Drought_Events_spring_start = [result_Flash_Drought_Events_spring_start, selected_many_FD_event_2D_spring_start];
    result_Flash_Drought_Events_spring_end = [result_Flash_Drought_Events_spring_end, selected_many_FD_event_2D_spring_end];
    % summer
    selected_many_FD_event_2D_summer = reshape(selected_many_FD_event_summer, 3900, col);
    selected_many_FD_event_2D_summer_start = reshape(selected_many_FD_event_summer_start, 3900, col);
    selected_many_FD_event_2D_summer_end = reshape(selected_many_FD_event_summer_end, 3900, col);

    result_Flash_Drought_Events_summer = [result_Flash_Drought_Events_summer, selected_many_FD_event_2D_summer];
    result_Flash_Drought_Events_summer_start = [result_Flash_Drought_Events_summer_start, selected_many_FD_event_2D_summer_start];
    result_Flash_Drought_Events_summer_end = [result_Flash_Drought_Events_summer_end, selected_many_FD_event_2D_summer_end];
    % autumn
    selected_many_FD_event_2D_autumn = reshape(selected_many_FD_event_autumn, 3900, col);
    selected_many_FD_event_2D_autumn_start = reshape(selected_many_FD_event_autumn_start, 3900, col);
    selected_many_FD_event_2D_autumn_end = reshape(selected_many_FD_event_autumn_end, 3900, col);

    result_Flash_Drought_Events_autumn = [result_Flash_Drought_Events_autumn, selected_many_FD_event_2D_autumn];
    result_Flash_Drought_Events_autumn_start = [result_Flash_Drought_Events_autumn_start, selected_many_FD_event_2D_autumn_start];
    result_Flash_Drought_Events_autumn_end = [result_Flash_Drought_Events_autumn_end, selected_many_FD_event_2D_autumn_end];
    % winter
    selected_many_FD_event_2D_winter = reshape(selected_many_FD_event_winter, 3900, col);
    selected_many_FD_event_2D_winter_start = reshape(selected_many_FD_event_winter_start, 3900, col);
    selected_many_FD_event_2D_winter_end = reshape(selected_many_FD_event_winter_end, 3900, col);

    result_Flash_Drought_Events_winter = [result_Flash_Drought_Events_winter, selected_many_FD_event_2D_winter];
    result_Flash_Drought_Events_winter_start = [result_Flash_Drought_Events_winter_start, selected_many_FD_event_2D_winter_start];
    result_Flash_Drought_Events_winter_end = [result_Flash_Drought_Events_winter_end, selected_many_FD_event_2D_winter_end];

    % 最慢的
    % spring
    selected_FD_event_point_slow_2D_spring = reshape(selected_FD_event_point_slow_spring, 3900, col);
    selected_FD_event_point_slow_2D_spring_start = reshape(selected_FD_event_point_slow_spring_start, 3900, col);
    selected_FD_event_point_slow_2D_spring_end = reshape(selected_FD_event_point_slow_spring_end, 3900, col);

    result_FD_Events_slow_spring = [result_FD_Events_slow_spring, selected_FD_event_point_slow_2D_spring];
    result_FD_Events_slow_spring_start = [result_FD_Events_slow_spring_start, selected_FD_event_point_slow_2D_spring_start];
    result_FD_Events_slow_spring_end = [result_FD_Events_slow_spring_end, selected_FD_event_point_slow_2D_spring_end];
    % summer
    selected_FD_event_point_slow_2D_summer = reshape(selected_FD_event_point_slow_summer, 3900, col);
    selected_FD_event_point_slow_2D_summer_start = reshape(selected_FD_event_point_slow_summer_start, 3900, col);
    selected_FD_event_point_slow_2D_summer_end = reshape(selected_FD_event_point_slow_summer_end, 3900, col);

    result_FD_Events_slow_summer = [result_FD_Events_slow_summer, selected_FD_event_point_slow_2D_summer];
    result_FD_Events_slow_summer_start = [result_FD_Events_slow_summer_start, selected_FD_event_point_slow_2D_summer_start];
    result_FD_Events_slow_summer_end = [result_FD_Events_slow_summer_end, selected_FD_event_point_slow_2D_summer_end];
    % autumn
    selected_FD_event_point_slow_2D_autumn = reshape(selected_FD_event_point_slow_autumn, 3900, col);
    selected_FD_event_point_slow_2D_autumn_start = reshape(selected_FD_event_point_slow_autumn_start, 3900, col);
    selected_FD_event_point_slow_2D_autumn_end = reshape(selected_FD_event_point_slow_autumn_end, 3900, col);

    result_FD_Events_slow_autumn = [result_FD_Events_slow_autumn, selected_FD_event_point_slow_2D_autumn];
    result_FD_Events_slow_autumn_start = [result_FD_Events_slow_autumn_start, selected_FD_event_point_slow_2D_autumn_start];
    result_FD_Events_slow_autumn_end = [result_FD_Events_slow_autumn_end, selected_FD_event_point_slow_2D_autumn_end];
    % winter
    selected_FD_event_point_slow_2D_winter = reshape(selected_FD_event_point_slow_winter, 3900, col);
    selected_FD_event_point_slow_2D_winter_start = reshape(selected_FD_event_point_slow_winter_start, 3900, col);
    selected_FD_event_point_slow_2D_winter_end = reshape(selected_FD_event_point_slow_winter_end, 3900, col);

    result_FD_Events_slow_winter = [result_FD_Events_slow_winter, selected_FD_event_point_slow_2D_winter];
    result_FD_Events_slow_winter_start = [result_FD_Events_slow_winter_start, selected_FD_event_point_slow_2D_winter_start];
    result_FD_Events_slow_winter_end = [result_FD_Events_slow_winter_end, selected_FD_event_point_slow_2D_winter_end];

    selected_SD_event_point_spring = zeros(3900*col,1);
    selected_SD_event_point_spring_start = zeros(3900*col,1);
    selected_SD_event_point_spring_end = zeros(3900*col,1);

    selected_SD_event_point_summer = zeros(3900*col,1);
    selected_SD_event_point_summer_start = zeros(3900*col,1);
    selected_SD_event_point_summer_end = zeros(3900*col,1);

    selected_SD_event_point_autumn = zeros(3900*col,1);
    selected_SD_event_point_autumn_start = zeros(3900*col,1);
    selected_SD_event_point_autumn_end = zeros(3900*col,1);

    selected_SD_event_point_winter = zeros(3900*col,1);
    selected_SD_event_point_winter_start = zeros(3900*col,1);
    selected_SD_event_point_winter_end = zeros(3900*col,1);

    selected_SD_event_point_slow_spring = zeros(3900*col,1);
    selected_SD_event_point_slow_spring_start = zeros(3900*col,1);
    selected_SD_event_point_slow_spring_end = zeros(3900*col,1);

    selected_SD_event_point_slow_summer = zeros(3900*col,1);
    selected_SD_event_point_slow_summer_start = zeros(3900*col,1);
    selected_SD_event_point_slow_summer_end = zeros(3900*col,1);

    selected_SD_event_point_slow_autumn = zeros(3900*col,1);
    selected_SD_event_point_slow_autumn_start = zeros(3900*col,1);
    selected_SD_event_point_slow_autumn_end = zeros(3900*col,1);

    selected_SD_event_point_slow_winter = zeros(3900*col,1);
    selected_SD_event_point_slow_winter_start = zeros(3900*col,1);
    selected_SD_event_point_slow_winter_end = zeros(3900*col,1);

    selected_many_SD_event_spring = cell(3900*col,1);
    selected_many_SD_event_spring_start = cell(3900*col,1);
    selected_many_SD_event_spring_end = cell(3900*col,1);

    selected_many_SD_event_summer = cell(3900*col,1);
    selected_many_SD_event_summer_start = cell(3900*col,1);
    selected_many_SD_event_summer_end = cell(3900*col,1);

    selected_many_SD_event_autumn = cell(3900*col,1);
    selected_many_SD_event_autumn_start = cell(3900*col,1);
    selected_many_SD_event_autumn_end = cell(3900*col,1);

    selected_many_SD_event_winter = cell(3900*col,1);
    selected_many_SD_event_winter_start = cell(3900*col,1);
    selected_many_SD_event_winter_end = cell(3900*col,1);

    for n = 1:numel(extracted_values_SD)
        
        FD_Events = extracted_values_FD{n};
        SD_Events = extracted_values_SD{n};

        if isempty(SD_Events)
            continue;
        end

        % 获取 NumObjects 字段的值
        SD_numObjects = SD_Events.NumObjects;

        if SD_numObjects == 0
            % 如果SD_numObjects为0，跳出循环
            continue;
        end

        DESI_contrast_values = extracted_rows_DESI(n, :);
        ESI_contrast_values = extracted_rows_ESI(n, :);

        position_index = positions(n);

        selected_SD_Event = [];  % 创建一个空的矩阵
        
        for e4 = 1:numel(SD_Events.PixelIdxList)
            SD_Event = SD_Events.PixelIdxList{e4};

            if isempty(SD_Event)
                continue
            end

            SD_Event_ESI_values = ESI_contrast_values(SD_Event);
            [SD_min_value, SD_min_index] = min(SD_Event_ESI_values);

            SD_Event_min_ESI_period = SD_Event(SD_min_index);
            SD_Event_first_ESI_period = SD_Event(1);

            if SD_Event_min_ESI_period <= 12 || SD_Event_min_ESI_period > 908
                % 最小值的周期在两端，则跳过当前循环
                continue;
            end

            SD_ESI_decrease_idx = SD_Event(1:SD_min_index);
            SD_DESI_decrease_idx = SD_ESI_decrease_idx(1:end-1);
            SD_Event_DESI_values = DESI_contrast_values(SD_DESI_decrease_idx);
            mean_SD_Event_DESI_values = mean(SD_Event_DESI_values);

            % 获取当前元素的第一个和最后一个值
            SD_Event_drought_start = SD_Event(1);
            SD_Event_drought_end = SD_Event(end);

            % 将这两个值存储为当前行
            current_row = [SD_Event_first_ESI_period, SD_Event_min_ESI_period, SD_min_value, mean_SD_Event_DESI_values, SD_Event_drought_start, SD_Event_drought_end];

            selected_SD_Event = [selected_SD_Event; current_row];
        end

        if isempty(selected_SD_Event)
            continue;
        end

        first_column2_spring = [];
        first_column2_spring_start = [];
        first_column2_spring_end = [];

        first_column2_summer = [];
        first_column2_summer_start = [];
        first_column2_summer_end = [];

        first_column2_autumn = [];
        first_column2_autumn_start = [];
        first_column2_autumn_end = [];

        first_column2_winter = [];
        first_column2_winter_start = [];
        first_column2_winter_end = [];

        % 记录添加的次数
        adding_count_spring = 0;
        adding_count_summer = 0;
        adding_count_autumn = 0;
        adding_count_winter = 0;

        % 按照第四列从小到大排列
        sorted_SD_Event = sortrows(selected_SD_Event, 4);
        
        % 循环遍历排序后的数组
        for s = 1:size(sorted_SD_Event, 1)
            % 获取第 i 行的数据
            drought_start = sorted_SD_Event(s, 1);
            drought_date = sorted_SD_Event(s, 2);

            % 获得前后 12 个索引号
            startIndex =  drought_date - 12;
            endIndex =  drought_date + 12;
            indices = startIndex:endIndex;

            overlapping_count = 0; % 记录重叠的次数

            for e5 = 1:numel(SD_Events.PixelIdxList)
                SD_Event = SD_Events.PixelIdxList{e5};
                
                % 检查SD_Event是否与indices有重叠
                if ~isempty(intersect(SD_Event, indices))
                    % 有重叠，增加重叠次数
                    overlapping_count = overlapping_count + 1;
                    
                    % 将该SD_Event添加到overlapping
                    only_one_overlapping = SD_Event;

                    % 如果重叠次数大于一个，直接进行下一次循环
                    if overlapping_count > 1
                        break;
                    end
                end
            end

            for e6 = 1:numel(FD_Events.PixelIdxList)

                FD_Event = FD_Events.PixelIdxList{e6};
                
                % 检查FD_Event是否与indices有重叠
                if ~isempty(intersect(FD_Event, indices))
                    % 有重叠，增加重叠次数
                    overlapping_count = overlapping_count + 1;

                    % 如果重叠次数大于一个，直接进行下一次循环
                    if overlapping_count > 1
                        break;
                    end
                end
            end

            % 如果重叠次数等于1
            if overlapping_count == 1
                % 检查drought_date是否包含在SD_Event中
                if ismember(drought_date, only_one_overlapping)
                    % 提取最小值对应的 FD_Event_min_ESI_period 值
                    result_SD_Event = drought_date;

                    % 获取开始的月份和年份
                    start_year = dateMatrix(drought_start, 1);
                    start_month = dateMatrix(drought_start, 2);
        
                    % 获取最严重时的月份和年份
                    end_year = dateMatrix(result_SD_Event, 1);
                    end_month = dateMatrix(result_SD_Event, 2);
                    
                    % 获取开始的季节和最严重时的季节
                    start_season = findSeason(start_month, interval_spring, interval_summer, interval_autumn, interval_winter);
                    end_season = findSeason(end_month, interval_spring, interval_summer, interval_autumn, interval_winter);

                    % 干旱开始
                    SD_drought_start = sorted_SD_Event(s, 5);
                    % 干旱结束
                    SD_drought_end = sorted_SD_Event(s, 6);
                    
                    % 根据不同季节执行相应操作
                    if isequal(start_season, end_season)
                        % 在相同的季节区间内
                        switch start_season
                            case 'spring'
                                if adding_count_spring == 0
                                    selected_SD_event_point_spring(position_index) = result_SD_Event;
                                    selected_SD_event_point_spring_start(position_index) = SD_drought_start;
                                    selected_SD_event_point_spring_end(position_index) = SD_drought_end;
                                    adding_count_spring = adding_count_spring + 1;
                                end
                                % 记录多个符合条件的事件
                                first_column2_spring(end+1) = result_SD_Event;
                                first_column2_spring_start(end+1) = SD_drought_start;
                                first_column2_spring_end(end+1) = SD_drought_end;
                            case 'summer'
                                if adding_count_summer == 0
                                    selected_SD_event_point_summer(position_index) = result_SD_Event;
                                    selected_SD_event_point_summer_start(position_index) = SD_drought_start;
                                    selected_SD_event_point_summer_end(position_index) = SD_drought_end;
                                    adding_count_summer = adding_count_summer + 1;
                                end
                                % 记录多个符合条件的事件
                                first_column2_summer(end+1) = result_SD_Event;
                                first_column2_summer_start(end+1) = SD_drought_start;
                                first_column2_summer_end(end+1) = SD_drought_end;
                            case 'autumn'
                                if adding_count_autumn == 0
                                    selected_SD_event_point_autumn(position_index) = result_SD_Event;
                                    selected_SD_event_point_autumn_start(position_index) = SD_drought_start;
                                    selected_SD_event_point_autumn_end(position_index) = SD_drought_end;
                                    adding_count_autumn = adding_count_autumn + 1;
                                end
                                % 记录多个符合条件的事件
                                first_column2_autumn(end+1) = result_SD_Event;
                                first_column2_autumn_start(end+1) = SD_drought_start;
                                first_column2_autumn_end(end+1) = SD_drought_end;
                            case 'winter'
                                if adding_count_winter == 0
                                    selected_SD_event_point_winter(position_index) = result_SD_Event;
                                    selected_SD_event_point_winter_start(position_index) = SD_drought_start;
                                    selected_SD_event_point_winter_end(position_index) = SD_drought_end;
                                    adding_count_winter = adding_count_winter + 1;
                                end
                                % 记录多个符合条件的事件
                                first_column2_winter(end+1) = result_SD_Event;
                                first_column2_winter_start(end+1) = SD_drought_start;
                                first_column2_winter_end(end+1) = SD_drought_end;
                            otherwise
                                disp('Unknown season.');
                                % 添加其他处理逻辑
                        end
                    else
                        % 不在相同的季节区间内
                        disp('Invalid season intervals for start_month and/or end_month.');
                        % 添加其他处理逻辑
                    end
                end
            end
        end

        if ~isempty(first_column2_spring)
            SD_slow_value_spring = first_column2_spring(end);
            SD_slow_value_spring_start = first_column2_spring_start(end);
            SD_slow_value_spring_end = first_column2_spring_end(end);

            selected_SD_event_point_slow_spring(position_index) = SD_slow_value_spring;
            selected_SD_event_point_slow_spring_start(position_index) = SD_slow_value_spring_start;
            selected_SD_event_point_slow_spring_end(position_index) = SD_slow_value_spring_end;

            selected_many_SD_event_spring{position_index} = first_column2_spring;
            selected_many_SD_event_spring_start{position_index} = first_column2_spring_start;
            selected_many_SD_event_spring_end{position_index} = first_column2_spring_end;
        end

        if ~isempty(first_column2_summer)
            SD_slow_value_summer = first_column2_summer(end);
            SD_slow_value_summer_start = first_column2_summer_start(end);
            SD_slow_value_summer_end = first_column2_summer_end(end);

            selected_SD_event_point_slow_summer(position_index) = SD_slow_value_summer;
            selected_SD_event_point_slow_summer_start(position_index) = SD_slow_value_summer_start;
            selected_SD_event_point_slow_summer_end(position_index) = SD_slow_value_summer_end;

            selected_many_SD_event_summer{position_index} = first_column2_summer;
            selected_many_SD_event_summer_start{position_index} = first_column2_summer_start;
            selected_many_SD_event_summer_end{position_index} = first_column2_summer_end;
        end

        if ~isempty(first_column2_autumn)
            SD_slow_value_autumn = first_column2_autumn(end);
            SD_slow_value_autumn_start = first_column2_autumn_start(end);
            SD_slow_value_autumn_end = first_column2_autumn_end(end);

            selected_SD_event_point_slow_autumn(position_index) = SD_slow_value_autumn;
            selected_SD_event_point_slow_autumn_start(position_index) = SD_slow_value_autumn_start;
            selected_SD_event_point_slow_autumn_end(position_index) = SD_slow_value_autumn_end;

            selected_many_SD_event_autumn{position_index} = first_column2_autumn;
            selected_many_SD_event_autumn_start{position_index} = first_column2_autumn_start;
            selected_many_SD_event_autumn_end{position_index} = first_column2_autumn_end;
        end

        if ~isempty(first_column2_winter)
            SD_slow_value_winter = first_column2_winter(end);
            SD_slow_value_winter_start = first_column2_winter_start(end);
            SD_slow_value_winter_end = first_column2_winter_end(end);

            selected_SD_event_point_slow_winter(position_index) = SD_slow_value_winter;
            selected_SD_event_point_slow_winter_start(position_index) = SD_slow_value_winter_start;
            selected_SD_event_point_slow_winter_end(position_index) = SD_slow_value_winter_end;

            selected_many_SD_event_winter{position_index} = first_column2_winter;
            selected_many_SD_event_winter_start{position_index} = first_column2_winter_start;
            selected_many_SD_event_winter_end{position_index} = first_column2_winter_end;
        end
    end

    % 最快的
    % spring
    selected_SD_event_point_2D_spring = reshape(selected_SD_event_point_spring, 3900, col);
    selected_SD_event_point_2D_spring_start = reshape(selected_SD_event_point_spring_start, 3900, col);
    selected_SD_event_point_2D_spring_end = reshape(selected_SD_event_point_spring_end, 3900, col);

    result_SD_Events_spring = [result_SD_Events_spring, selected_SD_event_point_2D_spring];
    result_SD_Events_spring_start = [result_SD_Events_spring_start, selected_SD_event_point_2D_spring_start];
    result_SD_Events_spring_end = [result_SD_Events_spring_end, selected_SD_event_point_2D_spring_end];
    % summer
    selected_SD_event_point_2D_summer = reshape(selected_SD_event_point_summer, 3900, col);
    selected_SD_event_point_2D_summer_start = reshape(selected_SD_event_point_summer_start, 3900, col);
    selected_SD_event_point_2D_summer_end = reshape(selected_SD_event_point_summer_end, 3900, col);

    result_SD_Events_summer = [result_SD_Events_summer, selected_SD_event_point_2D_summer];
    result_SD_Events_summer_start = [result_SD_Events_summer_start, selected_SD_event_point_2D_summer_start];
    result_SD_Events_summer_end = [result_SD_Events_summer_end, selected_SD_event_point_2D_summer_end];
    % autumn
    selected_SD_event_point_2D_autumn = reshape(selected_SD_event_point_autumn, 3900, col);
    selected_SD_event_point_2D_autumn_start = reshape(selected_SD_event_point_autumn_start, 3900, col);
    selected_SD_event_point_2D_autumn_end = reshape(selected_SD_event_point_autumn_end, 3900, col);

    result_SD_Events_autumn = [result_SD_Events_autumn, selected_SD_event_point_2D_autumn];
    result_SD_Events_autumn_start = [result_SD_Events_autumn_start, selected_SD_event_point_2D_autumn_start];
    result_SD_Events_autumn_end = [result_SD_Events_autumn_end, selected_SD_event_point_2D_autumn_end];
    % winter
    selected_SD_event_point_2D_winter = reshape(selected_SD_event_point_winter, 3900, col);
    selected_SD_event_point_2D_winter_start = reshape(selected_SD_event_point_winter_start, 3900, col);
    selected_SD_event_point_2D_winter_end = reshape(selected_SD_event_point_winter_end, 3900, col);

    result_SD_Events_winter = [result_SD_Events_winter, selected_SD_event_point_2D_winter];
    result_SD_Events_winter_start = [result_SD_Events_winter_start, selected_SD_event_point_2D_winter_start];
    result_SD_Events_winter_end = [result_SD_Events_winter_end, selected_SD_event_point_2D_winter_end];

    % 多个
    % spring
    selected_many_SD_event_2D_spring = reshape(selected_many_SD_event_spring, 3900, col);
    selected_many_SD_event_2D_spring_start = reshape(selected_many_SD_event_spring_start, 3900, col);
    selected_many_SD_event_2D_spring_end = reshape(selected_many_SD_event_spring_end, 3900, col);

    result_Slow_Drought_Events_spring = [result_Slow_Drought_Events_spring, selected_many_SD_event_2D_spring];
    result_Slow_Drought_Events_spring_start = [result_Slow_Drought_Events_spring_start, selected_many_SD_event_2D_spring_start];
    result_Slow_Drought_Events_spring_end = [result_Slow_Drought_Events_spring_end, selected_many_SD_event_2D_spring_end];
    % summer
    selected_many_SD_event_2D_summer = reshape(selected_many_SD_event_summer, 3900, col);
    selected_many_SD_event_2D_summer_start = reshape(selected_many_SD_event_summer_start, 3900, col);
    selected_many_SD_event_2D_summer_end = reshape(selected_many_SD_event_summer_end, 3900, col);

    result_Slow_Drought_Events_summer = [result_Slow_Drought_Events_summer, selected_many_SD_event_2D_summer];
    result_Slow_Drought_Events_summer_start = [result_Slow_Drought_Events_summer_start, selected_many_SD_event_2D_summer_start];
    result_Slow_Drought_Events_summer_end = [result_Slow_Drought_Events_summer_end, selected_many_SD_event_2D_summer_end];
    % autumn
    selected_many_SD_event_2D_autumn = reshape(selected_many_SD_event_autumn, 3900, col);
    selected_many_SD_event_2D_autumn_start = reshape(selected_many_SD_event_autumn_start, 3900, col);
    selected_many_SD_event_2D_autumn_end = reshape(selected_many_SD_event_autumn_end, 3900, col);

    result_Slow_Drought_Events_autumn = [result_Slow_Drought_Events_autumn, selected_many_SD_event_2D_autumn];
    result_Slow_Drought_Events_autumn_start = [result_Slow_Drought_Events_autumn_start, selected_many_SD_event_2D_autumn_start];
    result_Slow_Drought_Events_autumn_end = [result_Slow_Drought_Events_autumn_end, selected_many_SD_event_2D_autumn_end];
    % winter
    selected_many_SD_event_2D_winter = reshape(selected_many_SD_event_winter, 3900, col);
    selected_many_SD_event_2D_winter_start = reshape(selected_many_SD_event_winter_start, 3900, col);
    selected_many_SD_event_2D_winter_end = reshape(selected_many_SD_event_winter_end, 3900, col);

    result_Slow_Drought_Events_winter = [result_Slow_Drought_Events_winter, selected_many_SD_event_2D_winter];
    result_Slow_Drought_Events_winter_start = [result_Slow_Drought_Events_winter_start, selected_many_SD_event_2D_winter_start];
    result_Slow_Drought_Events_winter_end = [result_Slow_Drought_Events_winter_end, selected_many_SD_event_2D_winter_end];

    % 最慢的
    % spring
    selected_SD_event_point_slow_2D_spring = reshape(selected_SD_event_point_slow_spring, 3900, col);
    selected_SD_event_point_slow_2D_spring_start = reshape(selected_SD_event_point_slow_spring_start, 3900, col);
    selected_SD_event_point_slow_2D_spring_end = reshape(selected_SD_event_point_slow_spring_end, 3900, col);

    result_SD_Events_slow_spring = [result_SD_Events_slow_spring, selected_SD_event_point_slow_2D_spring];
    result_SD_Events_slow_spring_start = [result_SD_Events_slow_spring_start, selected_SD_event_point_slow_2D_spring_start];
    result_SD_Events_slow_spring_end = [result_SD_Events_slow_spring_end, selected_SD_event_point_slow_2D_spring_end];
    % summer
    selected_SD_event_point_slow_2D_summer = reshape(selected_SD_event_point_slow_summer, 3900, col);
    selected_SD_event_point_slow_2D_summer_start = reshape(selected_SD_event_point_slow_summer_start, 3900, col);
    selected_SD_event_point_slow_2D_summer_end = reshape(selected_SD_event_point_slow_summer_end, 3900, col);

    result_SD_Events_slow_summer = [result_SD_Events_slow_summer, selected_SD_event_point_slow_2D_summer];
    result_SD_Events_slow_summer_start = [result_SD_Events_slow_summer_start, selected_SD_event_point_slow_2D_summer_start];
    result_SD_Events_slow_summer_end = [result_SD_Events_slow_summer_end, selected_SD_event_point_slow_2D_summer_end];
    % autumn
    selected_SD_event_point_slow_2D_autumn = reshape(selected_SD_event_point_slow_autumn, 3900, col);
    selected_SD_event_point_slow_2D_autumn_start = reshape(selected_SD_event_point_slow_autumn_start, 3900, col);
    selected_SD_event_point_slow_2D_autumn_end = reshape(selected_SD_event_point_slow_autumn_end, 3900, col);

    result_SD_Events_slow_autumn = [result_SD_Events_slow_autumn, selected_SD_event_point_slow_2D_autumn];
    result_SD_Events_slow_autumn_start = [result_SD_Events_slow_autumn_start, selected_SD_event_point_slow_2D_autumn_start];
    result_SD_Events_slow_autumn_end = [result_SD_Events_slow_autumn_end, selected_SD_event_point_slow_2D_autumn_end];
    % winter
    selected_SD_event_point_slow_2D_winter = reshape(selected_SD_event_point_slow_winter, 3900, col);
    selected_SD_event_point_slow_2D_winter_start = reshape(selected_SD_event_point_slow_winter_start, 3900, col);
    selected_SD_event_point_slow_2D_winter_end = reshape(selected_SD_event_point_slow_winter_end, 3900, col);

    result_SD_Events_slow_winter = [result_SD_Events_slow_winter, selected_SD_event_point_slow_2D_winter];
    result_SD_Events_slow_winter_start = [result_SD_Events_slow_winter_start, selected_SD_event_point_slow_2D_winter_start];
    result_SD_Events_slow_winter_end = [result_SD_Events_slow_winter_end, selected_SD_event_point_slow_2D_winter_end];

    fprintf('已计算第 %d 块\n', t);
end

save_path_FD_spring = "E:\GZW\Drought_response\Season_N&P\Planted2\FD\spring";
save_path_FD_spring_periods = "E:\GZW\Drought_response\Season_N&P\Planted2\FD\spring\FD_spring_periods";

output_filename_FD_spring1 = "result_FD_Events_spring_8-10";
full_path_FD_spring1 = fullfile(save_path_FD_spring, output_filename_FD_spring1);
save(full_path_FD_spring1, 'result_FD_Events_spring', '-v7.3');

output_filename_FD_spring1_start = "result_FD_Events_spring_start_8-10";
full_path_FD_spring1_start = fullfile(save_path_FD_spring_periods, output_filename_FD_spring1_start);
save(full_path_FD_spring1_start, 'result_FD_Events_spring_start', '-v7.3');

output_filename_FD_spring1_end = "result_FD_Events_spring_end_8-10";
full_path_FD_spring1_end = fullfile(save_path_FD_spring_periods, output_filename_FD_spring1_end);
save(full_path_FD_spring1_end, 'result_FD_Events_spring_end', '-v7.3');

output_filename_FD_spring2 = "result_Flash_Drought_Events_spring_8-10";
full_path_FD_spring2 = fullfile(save_path_FD_spring, output_filename_FD_spring2);
save(full_path_FD_spring2, 'result_Flash_Drought_Events_spring', '-v7.3');

output_filename_FD_spring2_start = "result_Flash_Drought_Events_spring_start_8-10";
full_path_FD_spring2_start = fullfile(save_path_FD_spring_periods, output_filename_FD_spring2_start);
save(full_path_FD_spring2_start, 'result_Flash_Drought_Events_spring_start', '-v7.3');

output_filename_FD_spring2_end = "result_Flash_Drought_Events_spring_end_8-10";
full_path_FD_spring2_end = fullfile(save_path_FD_spring_periods, output_filename_FD_spring2_end);
save(full_path_FD_spring2_end, 'result_Flash_Drought_Events_spring_end', '-v7.3');

output_filename_FD_spring3 = "result_FD_Events_slow_spring_8-10";
full_path_FD_spring3 = fullfile(save_path_FD_spring, output_filename_FD_spring3);
save(full_path_FD_spring3, 'result_FD_Events_slow_spring', '-v7.3');

output_filename_FD_spring3_start = "result_FD_Events_slow_spring_start_8-10";
full_path_FD_spring3_start = fullfile(save_path_FD_spring_periods, output_filename_FD_spring3_start);
save(full_path_FD_spring3_start, 'result_FD_Events_slow_spring_start', '-v7.3');

output_filename_FD_spring3_end = "result_FD_Events_slow_spring_end_8-10";
full_path_FD_spring3_end = fullfile(save_path_FD_spring_periods, output_filename_FD_spring3_end);
save(full_path_FD_spring3_end, 'result_FD_Events_slow_spring_end', '-v7.3');

fprintf('已输出result_FD_Events_spring_8-10\n');

save_path_FD_summer = "E:\GZW\Drought_response\Season_N&P\Planted2\FD\summer";
save_path_FD_summer_periods = "E:\GZW\Drought_response\Season_N&P\Planted2\FD\summer\FD_summer_periods";

output_filename_FD_summer1 = "result_FD_Events_summer_8-10";
full_path_FD_summer1 = fullfile(save_path_FD_summer, output_filename_FD_summer1);
save(full_path_FD_summer1, 'result_FD_Events_summer', '-v7.3');

output_filename_FD_summer1_start = "result_FD_Events_summer_start_8-10";
full_path_FD_summer1_start = fullfile(save_path_FD_summer_periods, output_filename_FD_summer1_start);
save(full_path_FD_summer1_start, 'result_FD_Events_summer_start', '-v7.3');

output_filename_FD_summer1_end = "result_FD_Events_summer_end_8-10";
full_path_FD_summer1_end = fullfile(save_path_FD_summer_periods, output_filename_FD_summer1_end);
save(full_path_FD_summer1_end, 'result_FD_Events_summer_end', '-v7.3');

output_filename_FD_summer2 = "result_Flash_Drought_Events_summer_8-10";
full_path_FD_summer2 = fullfile(save_path_FD_summer, output_filename_FD_summer2);
save(full_path_FD_summer2, 'result_Flash_Drought_Events_summer', '-v7.3');

output_filename_FD_summer2_start = "result_Flash_Drought_Events_summer_start_8-10";
full_path_FD_summer2_start = fullfile(save_path_FD_summer_periods, output_filename_FD_summer2_start);
save(full_path_FD_summer2_start, 'result_Flash_Drought_Events_summer_start', '-v7.3');

output_filename_FD_summer2_end = "result_Flash_Drought_Events_summer_end_8-10";
full_path_FD_summer2_end = fullfile(save_path_FD_summer_periods, output_filename_FD_summer2_end);
save(full_path_FD_summer2_end, 'result_Flash_Drought_Events_summer_end', '-v7.3');

output_filename_FD_summer3 = "result_FD_Events_slow_summer_8-10";
full_path_FD_summer3 = fullfile(save_path_FD_summer, output_filename_FD_summer3);
save(full_path_FD_summer3, 'result_FD_Events_slow_summer', '-v7.3');

output_filename_FD_summer3_start = "result_FD_Events_slow_summer_start_8-10";
full_path_FD_summer3_start = fullfile(save_path_FD_summer_periods, output_filename_FD_summer3_start);
save(full_path_FD_summer3_start, 'result_FD_Events_slow_summer_start', '-v7.3');

output_filename_FD_summer3_end = "result_FD_Events_slow_summer_end_8-10";
full_path_FD_summer3_end = fullfile(save_path_FD_summer_periods, output_filename_FD_summer3_end);
save(full_path_FD_summer3_end, 'result_FD_Events_slow_summer_end', '-v7.3');

fprintf('已输出result_FD_Events_summer_8-10\n');

save_path_FD_autumn = "E:\GZW\Drought_response\Season_N&P\Planted2\FD\autumn";
save_path_FD_autumn_periods = "E:\GZW\Drought_response\Season_N&P\Planted2\FD\autumn\FD_autumn_periods";

output_filename_FD_autumn1 = "result_FD_Events_autumn_8-10";
full_path_FD_autumn1 = fullfile(save_path_FD_autumn, output_filename_FD_autumn1);
save(full_path_FD_autumn1, 'result_FD_Events_autumn', '-v7.3');

output_filename_FD_autumn1_start = "result_FD_Events_autumn_start_8-10";
full_path_FD_autumn1_start = fullfile(save_path_FD_autumn_periods, output_filename_FD_autumn1_start);
save(full_path_FD_autumn1_start, 'result_FD_Events_autumn_start', '-v7.3');

output_filename_FD_autumn1_end = "result_FD_Events_autumn_end_8-10";
full_path_FD_autumn1_end = fullfile(save_path_FD_autumn_periods, output_filename_FD_autumn1_end);
save(full_path_FD_autumn1_end, 'result_FD_Events_autumn_end', '-v7.3');

output_filename_FD_autumn2 = "result_Flash_Drought_Events_autumn_8-10";
full_path_FD_autumn2 = fullfile(save_path_FD_autumn, output_filename_FD_autumn2);
save(full_path_FD_autumn2, 'result_Flash_Drought_Events_autumn', '-v7.3');

output_filename_FD_autumn2_start = "result_Flash_Drought_Events_autumn_start_8-10";
full_path_FD_autumn2_start = fullfile(save_path_FD_autumn_periods, output_filename_FD_autumn2_start);
save(full_path_FD_autumn2_start, 'result_Flash_Drought_Events_autumn_start', '-v7.3');

output_filename_FD_autumn2_end = "result_Flash_Drought_Events_autumn_end_8-10";
full_path_FD_autumn2_end = fullfile(save_path_FD_autumn_periods, output_filename_FD_autumn2_end);
save(full_path_FD_autumn2_end, 'result_Flash_Drought_Events_autumn_end', '-v7.3');

output_filename_FD_autumn3 = "result_FD_Events_slow_autumn_8-10";
full_path_FD_autumn3 = fullfile(save_path_FD_autumn, output_filename_FD_autumn3);
save(full_path_FD_autumn3, 'result_FD_Events_slow_autumn', '-v7.3');

output_filename_FD_autumn3_start = "result_FD_Events_slow_autumn_start_8-10";
full_path_FD_autumn3_start = fullfile(save_path_FD_autumn_periods, output_filename_FD_autumn3_start);
save(full_path_FD_autumn3_start, 'result_FD_Events_slow_autumn_start', '-v7.3');

output_filename_FD_autumn3_end = "result_FD_Events_slow_autumn_end_8-10";
full_path_FD_autumn3_end = fullfile(save_path_FD_autumn_periods, output_filename_FD_autumn3_end);
save(full_path_FD_autumn3_end, 'result_FD_Events_slow_autumn_end', '-v7.3');

fprintf('已输出result_FD_Events_autumn_8-10\n');

save_path_FD_winter = "E:\GZW\Drought_response\Season_N&P\Planted2\FD\winter";
save_path_FD_winter_periods = "E:\GZW\Drought_response\Season_N&P\Planted2\FD\winter\FD_winter_periods";

output_filename_FD_winter1 = "result_FD_Events_winter_8-10";
full_path_FD_winter1 = fullfile(save_path_FD_winter, output_filename_FD_winter1);
save(full_path_FD_winter1, 'result_FD_Events_winter', '-v7.3');

output_filename_FD_winter1_start = "result_FD_Events_winter_start_8-10";
full_path_FD_winter1_start = fullfile(save_path_FD_winter_periods, output_filename_FD_winter1_start);
save(full_path_FD_winter1_start, 'result_FD_Events_winter_start', '-v7.3');

output_filename_FD_winter1_end = "result_FD_Events_winter_end_8-10";
full_path_FD_winter1_end = fullfile(save_path_FD_winter_periods, output_filename_FD_winter1_end);
save(full_path_FD_winter1_end, 'result_FD_Events_winter_end', '-v7.3');

output_filename_FD_winter2 = "result_Flash_Drought_Events_winter_8-10";
full_path_FD_winter2 = fullfile(save_path_FD_winter, output_filename_FD_winter2);
save(full_path_FD_winter2, 'result_Flash_Drought_Events_winter', '-v7.3');

output_filename_FD_winter2_start = "result_Flash_Drought_Events_winter_start_8-10";
full_path_FD_winter2_start = fullfile(save_path_FD_winter_periods, output_filename_FD_winter2_start);
save(full_path_FD_winter2_start, 'result_Flash_Drought_Events_winter_start', '-v7.3');

output_filename_FD_winter2_end = "result_Flash_Drought_Events_winter_end_8-10";
full_path_FD_winter2_end = fullfile(save_path_FD_winter_periods, output_filename_FD_winter2_end);
save(full_path_FD_winter2_end, 'result_Flash_Drought_Events_winter_end', '-v7.3');

output_filename_FD_winter3 = "result_FD_Events_slow_winter_8-10";
full_path_FD_winter3 = fullfile(save_path_FD_winter, output_filename_FD_winter3);
save(full_path_FD_winter3, 'result_FD_Events_slow_winter', '-v7.3');

output_filename_FD_winter3_start = "result_FD_Events_slow_winter_start_8-10";
full_path_FD_winter3_start = fullfile(save_path_FD_winter_periods, output_filename_FD_winter3_start);
save(full_path_FD_winter3_start, 'result_FD_Events_slow_winter_start', '-v7.3');

output_filename_FD_winter3_end = "result_FD_Events_slow_winter_end_8-10";
full_path_FD_winter3_end = fullfile(save_path_FD_winter_periods, output_filename_FD_winter3_end);
save(full_path_FD_winter3_end, 'result_FD_Events_slow_winter_end', '-v7.3');

fprintf('已输出result_FD_Events_winter_8-10\n');

% SD
save_path_SD_spring = "E:\GZW\Drought_response\Season_N&P\Planted2\SD\spring";
save_path_SD_spring_periods = "E:\GZW\Drought_response\Season_N&P\Planted2\SD\spring\SD_spring_periods";

output_filename_SD_spring1 = "result_SD_Events_spring_8-10";
full_path_SD_spring1 = fullfile(save_path_SD_spring, output_filename_SD_spring1);
save(full_path_SD_spring1, 'result_SD_Events_spring', '-v7.3');

output_filename_SD_spring1_start = "result_SD_Events_spring_start_8-10";
full_path_SD_spring1_start = fullfile(save_path_SD_spring_periods, output_filename_SD_spring1_start);
save(full_path_SD_spring1_start, 'result_SD_Events_spring_start', '-v7.3');

output_filename_SD_spring1_end = "result_SD_Events_spring_end_8-10";
full_path_SD_spring1_end = fullfile(save_path_SD_spring_periods, output_filename_SD_spring1_end);
save(full_path_SD_spring1_end, 'result_SD_Events_spring_end', '-v7.3');

output_filename_SD_spring2 = "result_Slow_Drought_Events_spring_8-10";
full_path_SD_spring2 = fullfile(save_path_SD_spring, output_filename_SD_spring2);
save(full_path_SD_spring2, 'result_Slow_Drought_Events_spring', '-v7.3');

output_filename_SD_spring2_start = "result_Slow_Drought_Events_spring_start_8-10";
full_path_SD_spring2_start = fullfile(save_path_SD_spring_periods, output_filename_SD_spring2_start);
save(full_path_SD_spring2_start, 'result_Slow_Drought_Events_spring_start', '-v7.3');

output_filename_SD_spring2_end = "result_Slow_Drought_Events_spring_end_8-10";
full_path_SD_spring2_end = fullfile(save_path_SD_spring_periods, output_filename_SD_spring2_end);
save(full_path_SD_spring2_end, 'result_Slow_Drought_Events_spring_end', '-v7.3');

output_filename_SD_spring3 = "result_SD_Events_slow_spring_8-10";
full_path_SD_spring3 = fullfile(save_path_SD_spring, output_filename_SD_spring3);
save(full_path_SD_spring3, 'result_SD_Events_slow_spring', '-v7.3');

output_filename_SD_spring3_start = "result_SD_Events_slow_spring_start_8-10";
full_path_SD_spring3_start = fullfile(save_path_SD_spring_periods, output_filename_SD_spring3_start);
save(full_path_SD_spring3_start, 'result_SD_Events_slow_spring_start', '-v7.3');

output_filename_SD_spring3_end = "result_SD_Events_slow_spring_end_8-10";
full_path_SD_spring3_end = fullfile(save_path_SD_spring_periods, output_filename_SD_spring3_end);
save(full_path_SD_spring3_end, 'result_SD_Events_slow_spring_end', '-v7.3');

fprintf('已输出result_SD_Events_spring_8-10\n');

save_path_SD_summer = "E:\GZW\Drought_response\Season_N&P\Planted2\SD\summer";
save_path_SD_summer_periods = "E:\GZW\Drought_response\Season_N&P\Planted2\SD\summer\SD_summer_periods";

output_filename_SD_summer1 = "result_SD_Events_summer_8-10";
full_path_SD_summer1 = fullfile(save_path_SD_summer, output_filename_SD_summer1);
save(full_path_SD_summer1, 'result_SD_Events_summer', '-v7.3');

output_filename_SD_summer1_start = "result_SD_Events_summer_start_8-10";
full_path_SD_summer1_start = fullfile(save_path_SD_summer_periods, output_filename_SD_summer1_start);
save(full_path_SD_summer1_start, 'result_SD_Events_summer_start', '-v7.3');

output_filename_SD_summer1_end = "result_SD_Events_summer_end_8-10";
full_path_SD_summer1_end = fullfile(save_path_SD_summer_periods, output_filename_SD_summer1_end);
save(full_path_SD_summer1_end, 'result_SD_Events_summer_end', '-v7.3');

output_filename_SD_summer2 = "result_Slow_Drought_Events_summer_8-10";
full_path_SD_summer2 = fullfile(save_path_SD_summer, output_filename_SD_summer2);
save(full_path_SD_summer2, 'result_Slow_Drought_Events_summer', '-v7.3');

output_filename_SD_summer2_start = "result_Slow_Drought_Events_summer_start_8-10";
full_path_SD_summer2_start = fullfile(save_path_SD_summer_periods, output_filename_SD_summer2_start);
save(full_path_SD_summer2_start, 'result_Slow_Drought_Events_summer_start', '-v7.3');

output_filename_SD_summer2_end = "result_Slow_Drought_Events_summer_end_8-10";
full_path_SD_summer2_end = fullfile(save_path_SD_summer_periods, output_filename_SD_summer2_end);
save(full_path_SD_summer2_end, 'result_Slow_Drought_Events_summer_end', '-v7.3');

output_filename_SD_summer3 = "result_SD_Events_slow_summer_8-10";
full_path_SD_summer3 = fullfile(save_path_SD_summer, output_filename_SD_summer3);
save(full_path_SD_summer3, 'result_SD_Events_slow_summer', '-v7.3');

output_filename_SD_summer3_start = "result_SD_Events_slow_summer_start_8-10";
full_path_SD_summer3_start = fullfile(save_path_SD_summer_periods, output_filename_SD_summer3_start);
save(full_path_SD_summer3_start, 'result_SD_Events_slow_summer_start', '-v7.3');

output_filename_SD_summer3_end = "result_SD_Events_slow_summer_end_8-10";
full_path_SD_summer3_end = fullfile(save_path_SD_summer_periods, output_filename_SD_summer3_end);
save(full_path_SD_summer3_end, 'result_SD_Events_slow_summer_end', '-v7.3');

fprintf('已输出result_SD_Events_summer_8-10\n');

save_path_SD_autumn = "E:\GZW\Drought_response\Season_N&P\Planted2\SD\autumn";
save_path_SD_autumn_periods = "E:\GZW\Drought_response\Season_N&P\Planted2\SD\autumn\SD_autumn_periods";

output_filename_SD_autumn1 = "result_SD_Events_autumn_8-10";
full_path_SD_autumn1 = fullfile(save_path_SD_autumn, output_filename_SD_autumn1);
save(full_path_SD_autumn1, 'result_SD_Events_autumn', '-v7.3');

output_filename_SD_autumn1_start = "result_SD_Events_autumn_start_8-10";
full_path_SD_autumn1_start = fullfile(save_path_SD_autumn_periods, output_filename_SD_autumn1_start);
save(full_path_SD_autumn1_start, 'result_SD_Events_autumn_start', '-v7.3');

output_filename_SD_autumn1_end = "result_SD_Events_autumn_end_8-10";
full_path_SD_autumn1_end = fullfile(save_path_SD_autumn_periods, output_filename_SD_autumn1_end);
save(full_path_SD_autumn1_end, 'result_SD_Events_autumn_end', '-v7.3');

output_filename_SD_autumn2 = "result_Slow_Drought_Events_autumn_8-10";
full_path_SD_autumn2 = fullfile(save_path_SD_autumn, output_filename_SD_autumn2);
save(full_path_SD_autumn2, 'result_Slow_Drought_Events_autumn', '-v7.3');

output_filename_SD_autumn2_start = "result_Slow_Drought_Events_autumn_start_8-10";
full_path_SD_autumn2_start = fullfile(save_path_SD_autumn_periods, output_filename_SD_autumn2_start);
save(full_path_SD_autumn2_start, 'result_Slow_Drought_Events_autumn_start', '-v7.3');

output_filename_SD_autumn2_end = "result_Slow_Drought_Events_autumn_end_8-10";
full_path_SD_autumn2_end = fullfile(save_path_SD_autumn_periods, output_filename_SD_autumn2_end);
save(full_path_SD_autumn2_end, 'result_Slow_Drought_Events_autumn_end', '-v7.3');

output_filename_SD_autumn3 = "result_SD_Events_slow_autumn_8-10";
full_path_SD_autumn3 = fullfile(save_path_SD_autumn, output_filename_SD_autumn3);
save(full_path_SD_autumn3, 'result_SD_Events_slow_autumn', '-v7.3');

output_filename_SD_autumn3_start = "result_SD_Events_slow_autumn_start_8-10";
full_path_SD_autumn3_start = fullfile(save_path_SD_autumn_periods, output_filename_SD_autumn3_start);
save(full_path_SD_autumn3_start, 'result_SD_Events_slow_autumn_start', '-v7.3');

output_filename_SD_autumn3_end = "result_SD_Events_slow_autumn_end_8-10";
full_path_SD_autumn3_end = fullfile(save_path_SD_autumn_periods, output_filename_SD_autumn3_end);
save(full_path_SD_autumn3_end, 'result_SD_Events_slow_autumn_end', '-v7.3');

fprintf('已输出result_SD_Events_autumn_8-10\n');

save_path_SD_winter = "E:\GZW\Drought_response\Season_N&P\Planted2\SD\winter";
save_path_SD_winter_periods = "E:\GZW\Drought_response\Season_N&P\Planted2\SD\winter\SD_winter_periods";

output_filename_SD_winter1 = "result_SD_Events_winter_8-10";
full_path_SD_winter1 = fullfile(save_path_SD_winter, output_filename_SD_winter1);
save(full_path_SD_winter1, 'result_SD_Events_winter', '-v7.3');

output_filename_SD_winter1_start = "result_SD_Events_winter_start_8-10";
full_path_SD_winter1_start = fullfile(save_path_SD_winter_periods, output_filename_SD_winter1_start);
save(full_path_SD_winter1_start, 'result_SD_Events_winter_start', '-v7.3');

output_filename_SD_winter1_end = "result_SD_Events_winter_end_8-10";
full_path_SD_winter1_end = fullfile(save_path_SD_winter_periods, output_filename_SD_winter1_end);
save(full_path_SD_winter1_end, 'result_SD_Events_winter_end', '-v7.3');

output_filename_SD_winter2 = "result_Slow_Drought_Events_winter_8-10";
full_path_SD_winter2 = fullfile(save_path_SD_winter, output_filename_SD_winter2);
save(full_path_SD_winter2, 'result_Slow_Drought_Events_winter', '-v7.3');

output_filename_SD_winter2_start = "result_Slow_Drought_Events_winter_start_8-10";
full_path_SD_winter2_start = fullfile(save_path_SD_winter_periods, output_filename_SD_winter2_start);
save(full_path_SD_winter2_start, 'result_Slow_Drought_Events_winter_start', '-v7.3');

output_filename_SD_winter2_end = "result_Slow_Drought_Events_winter_end_8-10";
full_path_SD_winter2_end = fullfile(save_path_SD_winter_periods, output_filename_SD_winter2_end);
save(full_path_SD_winter2_end, 'result_Slow_Drought_Events_winter_end', '-v7.3');

output_filename_SD_winter3 = "result_SD_Events_slow_winter_8-10";
full_path_SD_winter3 = fullfile(save_path_SD_winter, output_filename_SD_winter3);
save(full_path_SD_winter3, 'result_SD_Events_slow_winter', '-v7.3');

output_filename_SD_winter3_start = "result_SD_Events_slow_winter_start_8-10";
full_path_SD_winter3_start = fullfile(save_path_SD_winter_periods, output_filename_SD_winter3_start);
save(full_path_SD_winter3_start, 'result_SD_Events_slow_winter_start', '-v7.3');

output_filename_SD_winter3_end = "result_SD_Events_slow_winter_end_8-10";
full_path_SD_winter3_end = fullfile(save_path_SD_winter_periods, output_filename_SD_winter3_end);
save(full_path_SD_winter3_end, 'result_SD_Events_slow_winter_end', '-v7.3');

fprintf('已输出result_SD_Events_winter_8-10\n');


% 季节判断函数
function season = findSeason(month, interval_spring, interval_summer, interval_autumn, interval_winter)
    if ismember(month, interval_spring)
        season = 'spring';
    elseif ismember(month, interval_summer)
        season = 'summer';
    elseif ismember(month, interval_autumn)
        season = 'autumn';
    elseif ismember(month, interval_winter)
        season = 'winter';
    else
        season = 'unknown';
    end
end