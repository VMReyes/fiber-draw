function [X, Y] = load_good_fiber_data(strDataPath, x_columns, y_columns, loBFD, hiBFD)

all_files = dir(strDataPath);
headers = readcell(strcat(strDataPath,all_files(3).name),'Range','1:1');
for curr_file_ind = 3:16
    x_data = [];
    y_data = [];
    for y_column_name = y_columns
        y_column_letter = xlsColNum2Str(find(strcmpi(headers, y_column_name)==1));
        y_data{end+1} = readmatrix(strcat(strDataPath,all_files(curr).name), 'Range', sprintf("%s:%s", y_column_letter, y_column_letter));
    end
    for x_column_name = x_columns
        x_column_letter = xlsColNum2Str(find(strcmpi(headers, x_column_name)==1));
        x_data{end+1} = readmatrix([strDataPath all_files(curr).name], 'Range', sprintf("%s:%s", x_column_letter, x_column_letter));
    end

    bfd = readmatrix([strDataPath all_files(curr_file_ind).name], 'Range', 'C:C');
    good_fiber = readmatrix([strDataPath all_files(curr_file_ind).name], 'Range', 'Q:Q');

    x_subbatches = [];
    y_subbatches = [];
    % want continuous data within range and on rising edge of good_fiber
    prev_good_fiber = 0;
    for i = range(length(bfd))
        x_subbatch = [];
        y_subbatch = [];
        if (prev_good_fiber == 0) && (good_fiber(i) == 1)
            x_subbatch{end+1} = x_data(i,:);
            y_subbatch{end+1} = y_data(i,:);
        end
        if (prev_good_fiber == 1) && (good_fiber(i) == 1)
            x_subbatch{end+1} = x_data(i,:);
            y_subbatch{end+1} = y_data(i,:);
        end
        if (prev_good_fiber == 1) && (good_fiber == 0)
            x_subbatches{end+1} = x_subbatch;
            y_subbatches{end+1} = y_subbatch;
            x_subbatch = [];
            y_subbatch = [];
        end
        prev_good_fiber = good_fiber(i);
    end
end