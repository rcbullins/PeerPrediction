function [dev_final, devControl_final] = PP_Load_and_Concatenate(data_path)
%Purpose: This function will load all data files from the PP function.
%         And all files will be concatenated.
%Inputs: data path

%Outputs: dev_final : may be in one matrix if concatenated, or may be in
%                     array if not concatenated
%         devControl_final

% Created: 3/23/20 by Reagan Bullins

%% go to data folder
    cd(data_path)

%% load in data & concatenat if prompted
%find folder name
    [~,folder_name,~] = fileparts(data_path);
%find directories for all files with folder name
    files_struct = what(folder_name)
%find all mat files in main folder name
    mat_files = files_struct(1).mat;
%concat all file names into matrix
    mat_files = cat(1,mat_files{:});  
    [num_of_files,~] = size(mat_files);
    
    resultDev = {};
    resultDevControl = {};
%for every mat file, load the mat file, and save it to an array
    for iResult = 1:num_of_files
        load(mat_files(iResult,:))
        resultDev{iResult} = dev
        resultDevControl{iResult} = devControl
        clear dev
        clear devControl
    end
%concat all mat file arrays together
        dev_final = cat(2, resultDev{:});
        devControl_final = cat(2, resultDevControl{:});

end
  
