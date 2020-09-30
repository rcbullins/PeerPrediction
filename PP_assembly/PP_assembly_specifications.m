%% PP_assembly_specifications

%Add Basepath for all code and data used
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
%Define Recording Session Name
    %session_name = 'm115_191203_152410_n';
    session_name = 'u19_200313_155505'; %NUM 1 -running
    %session_name = 'u19_200310_135409'; %NUM 2
    %session_name = 'u19_200313_120452'; %NUM 3
    %session_name = 'u21_200305_153604'; %NUM 4 -runninng
    %session_name = 'u21_200309_142534'; %NUM 5
    %session_name = 'u26_200306_172032'; %NUM 6
%Deine DataPath that contains list of session names;
    data_path = [basepath 'PP_RSC_Data\' session_name];
%Add Paths
    addpath(genpath(data_path));
    addpath(genpath([basepath 'Code\']));
    addpath(genpath([basepath 'Sam_Code\']));
    addpath(genpath([basepath 'buzcode-dev\']));
  
%% Get epochs of running vs no running
    cd(data_path)
    % speed of the animal to be over 1.5cm/s for RUN
    
    % define parameters - length of track, etc
        rad_disk = 26; %cm
        circum_disk = 163.3628; %(2*pi*r) cm
        load([session_name '_analogin.mat']);
    % downsample 30000ms to 1 second (oh yes)
        ds_analogin_pos = analogin.pos(1:30000:length(analogin.pos)); %analogin once every second
        
                %plot for reference to see downsampling effect
            %          plot(1:length(analogin.pos), analogin.pos, 'k');
            %          hold on;
            %          plot(1:30:length(analogin.pos), ds_analogin_pos(1,:));

    % smooth downsampled data - averaging
        smooth_ds_pos = smoothdata(ds_analogin_pos);
                %plot for reference to see smoothing effect
%                     plot(1:length(analogin.pos), analogin.pos, 'k');
%                     hold on;
%                     plot(1:30:length(analogin.pos), smooth_ds_pos(1,:));
%                 
        smooth_ds_pos = round(smooth_ds_pos,3); %round to 3rd decimal
   
    % find the min and max voltage
        max_volt = max(smooth_ds_pos);
        min_volt = min(smooth_ds_pos);
    % make .001 intervals between min and max voltage
        pos_convert = (min_volt:.001:max_volt);
        pos_convert = round(pos_convert,3); %VERY important 
    % want to assign each distance point to voltage point
        convert_ratio =  circum_disk/length(pos_convert);
        pos_convert(1,:) = pos_convert; %first row equals voltage, second row equals position
        pos_convert(2,:) = (0.001:convert_ratio:circum_disk); %think .001 works?... idk

    % whenever row 1 equals this, make it same index row 2
    % initiate a vector: this will be the newly created position vector per 1 ms
        smooth_ds_pos_transition = NaN(1, length(smooth_ds_pos)); % new vector length of volt vector
        %for every possible voltage point, see if there is coinciding
        %position points
        for ivolt = 1:length(pos_convert(1,:))
            idx_change_volt = find(smooth_ds_pos == pos_convert(1, ivolt)); %find idx where voltage equals the specified voltage
            smooth_ds_pos_transition(1, idx_change_volt) =  pos_convert(2, ivolt); %convert these values to pos
        end
            smooth_ds_pos_convert = smooth_ds_pos_transition;
        diff_filtered_pos = abs(diff(smooth_ds_pos_convert));
    % calculate velocity 
        %units depends on downsampled (if to 1000 per second, this is
        %cm/ms)
        vel_pos = diff_filtered_pos; %cm diff per second
    % smooth velocity
  % ________________________________________________________________
        smth_vel_pos = smoothdata(vel_pos, 'sgolay', 60) 
        smth_vel1 = smoothdata(vel_pos, 'movmean',10)
        smth_vel2 = smoothdata(vel_pos, 'movmean',30)
        smth_vel3 = smoothdata(vel_pos,'movmean',60)
  %__________________________________________________________________
    % find velocity above arbitrary moving threshold
        above_move_thres = find(vel_pos(1,:) >= 1.5); 
   
        %see how long in total there is movement
        move_tot_sec = length(above_move_thres)/1000 %i think?? need this divide
        plot(length(vel_pos(1,1000)), vel_pos(1,1000));

        histogram(vel_pos)
 %% Find Epochs of running vs not running
      
    
    