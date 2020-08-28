%% PP_assembly_specifications

%Add Basepath for all code and data used
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
%Define Recording Session Name
    %session_name = 'm115_191203_152410_n';
    session_name = 'u19_200313_155505';
%Deine DataPath that contains list of session names;
    data_path = [basepath 'PP_RSC_Data\' session_name];
%Add Paths
    addpath(genpath(data_path));
    addpath(genpath([basepath 'Code\']));
    addpath(genpath([basepath 'Sam_Code\']));
    addpath(genpath([basepath 'buzcode-dev\']));

%% Get epochs of running vs no running
    % speed of the animal to be over 1.5cm/s for RUN
    rad_disk = 26 %cm
    circum_disk = 163.3628 %(2*pi*r) cm
    load([session_name '_analogin.mat']);
    % downsample 30000 to 1000 :)
    ds_analogin_pos = analogin.pos(1:30:length(analogin.pos));
    
    %plot for reference to see downsampling effect
    %     plot(1:length(analogin.pos), analogin.pos, 'k');
    %     hold on;
    %     plot(1:30:length(analogin.pos), ds_analogin_pos(1,:));

    % smooth downsampled data - averaging
    smooth_ds_pos = smoothdata(ds_analogin_pos);
    %plot for reference to see smoothing effect
%         plot(1:length(analogin.pos), analogin.pos, 'k');
%         hold on;
%         plot(1:30:length(analogin.pos), smooth_ds_pos(1,:));
%     
    % find the difference in distance between consecutive points
    % first, convert the voltage to distance points
    max_volt = max(smooth_ds_pos);
    min_volt = min(smooth_ds_pos);
    
    smooth_ds_pos = round(smooth_ds_pos,3);
    pos_convert = unique(smooth_ds_pos);
    
%     volt_diff = max_volt - min_volt;
%     interval_size = volt_diff/circum_disk;

    convert_ratio =  circum_disk/length(pos_convert);
    pos_convert(1,:) = pos_convert;
    pos_convert(2,:) = (0.001:convert_ratio:circum_disk); %think .001 works?... idk
    
    %WORK HERE::::need to say, whenever row 1 equals this, make it same index row 2
    smooth_ds_pos_convert =  
    
    diff_filtered_pos = diff(smooth_ds_pos_convert);
    vel_pos = diff_filtered_pos/.001; %cm/ms
    above_move_thres = find(vel_pos(1,:) >= 1.5); %seemes to be too high
    plot(length(vel_pos(1,:)), vel_pos(1,:));
    
    
    
    