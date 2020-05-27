function [dev_final, devControl_final] = PP_Load_and_Concatenate(data_path)
%Purpose: This function will load all OR specific result data files from the PP function
%         A prompt will be given to choose to load all data or specify. Additionally, the code
%         will prompt if you would like to concatenate the data or not.

%Inputs: data path

%Outputs: dev_final : may be in one matrix if concatenated, or may be in
%                     array if not concatenated
%         devControl_final

% Created: 3/23/20 by Reagan Bullins

%% go to data folder
cd(data_path)

%% load in data & concatenat if prompted 

prompt = 'Run all data in folder or specify? All OR Specify:'
choice = input(prompt,'s')

if strcmp(choice, 'All')
    prompt2 = 'Concatenate all data? Yes OR No:'
    concat_choice = input(prompt2,'s');
    
    [~,folder_name,~] = fileparts(data_path);
    files_struct = what(folder_name)
    mat_files = files_struct(1).mat;
    mat_files = cat(1,mat_files{:});  
    [num_of_files,~] = size(mat_files);
    
    resultDev = {};
    resultDevControl = {};

    for iResult = 1:num_of_files
        load(mat_files(iResult,:))
        resultDev{iResult} = dev
        resultDevControl{iResult} = devControl
        clear dev
        clear devControl
    end
    
    if strcmp(concat_choice, 'Yes')
        dev_final = cat(2, resultDev{:});
        devControl_final = cat(2, resultDevControl{:});
    
    else 
        dev_final = resultDev;
        devControl_final = resultDevControl;
    end
        
        
else if strcmp(choice, 'Specify')
    prompt2 = 'How many files do you want to load?:'
    num_files_to_load = input(prompt2);
    
    if num_files_to_load > 1
        prompt3 = 'Would you like them concatenated together? Yes OR No:'
        concat_choice = input(prompt3, 's')
    end
    
    chosen_files = {}
    for ifile = 1:num_files_to_load
        prompt_file = ['List File ' num2str(ifile) 'with ext .mat you want to run:']
        chosen_files{ifile} = input(prompt_file,'s')
    end
    
    %chosen_files = cat(1,chosen_files{:}); 
    resultDev = {};
    resultDevControl = {};

    for iResult = 1:length(chosen_files)
        load(chosen_files{iResult})
        resultDev{iResult} = dev
        resultDevControl{iResult} = devControl
        clear dev
        clear devControl
    end
    
    if strcmp(concat_choice, 'Yes')
       dev_final = cat(2, resultDev{:});
       devControl_final = cat(2, resultDevControl{:});
    else
        dev_final = resultDev;
        devControl_final = resultDevControl;
    end
    

end 
end
  
%%
% resultDev(zeros(151, length(mat_files)*num_cell_pairs)); %151 is timewindows, default 151 given by output of PP
% resultDevControl(zeros(151, length(mat_files)*num_cell_pairs), num_control_trials);
%%%
% load(mat_files(iResult),'-mat')
% resultDev(iResult) = dev
% resultDevControl(iResult) = devControl
% clear dev
% clear devControl
% %%%
% result1 = load(results_a, '-mat')
% dev1 = dev
% devControl1 = devControl
% clear dev
% clear devControl
% 
% result2 = load(results_b, '-mat')
% dev2 = dev
% devControl2 = devControl
% clear dev
% clear devControl
%% Concatenate mat files
% dev_cat =(1, resultDev{:});
% devControl_cat = (1, resultDevControl{:});

% dev_cat = [dev1,dev2]
% devControl_cat = [devControl1, devControl2]


