function [ratio_strength, dev_min, iweak_pair, istrong_pair] = PP_AssemblyStrength (dev, devControl, data_folder)
%Purpose: To find assembly strength of each cell pair, see how well each
%         pair is more or less correlated to one another

%Inputs: Ouputs of PP function including dev and devControl

%Outputs: ratio_Strength (assembly strength per cell)
%         dev_min (minimum deviance value for each pair)

%Created: 3/30/20 by Reagan Bullins

[~,number_of_pairs] = size(dev);
dev_min = zeros(1, number_of_pairs);
devControl_min = zeros(1, number_of_pairs);
ratio_strength = zeros(1, number_of_pairs);

for ipair = 1:number_of_pairs
    dev_min(ipair) = min(dev(:,ipair));
    devControl_min(ipair) = min(mean(devControl(:,ipair,:),3)); % or mean across 3D--> mean(devControl,3);
    ratio_strength(ipair) = (dev_min(ipair) - mean(dev(:,ipair))) ./ ...
                            (devControl_min(ipair) - mean(mean(devControl(:,ipair,:),3)));
end

if strcmp(data_folder, 'pp_batch')
    iweak_pair = find(ratio_strength < 3.50)
    istrong_pair = find(ratio_strength >= 3.50)
elseif strcmp(data_folder, 'pp_poisson')
    iweak_pair = find(ratio_strength < 350)
    istrong_pair = find(ratio_strength >= 350)
end

%weak_strength = ratio_strength(iweak_pair)
%strong_strength = ratio_strength(istrong_pair)

end