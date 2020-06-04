function [dev_smoothed, ratio_strength, dev_min, iweak_pair, istrong_pair] = PP_AssemblyStrength (dev, devControl, data_folder, choice_smooth)
%Purpose: To find assembly strength of each cell pair, see how well each
%         pair is more or less correlated to one another

%Inputs: Ouputs of PP function including dev and devControl

%Outputs: ratio_Strength & Smoothed(assembly strength per cell)
%         dev_min & Smoothed (minimum deviance value for each pair)
%         weak_pairs (assemb strength)
%         strong_pairs (assemb strength)
%         dev_smoothed

%Created: 03/30/20 by Reagan Bullins
%Updated: 06/03/20 by Reagan Bullins
  [~,number_of_pairs] = size(dev);
if choice_smooth == 0
    dev_min = zeros(1, number_of_pairs);
    ratio_strength = zeros(1, number_of_pairs);
    devControl_min = zeros(1, number_of_pairs);
    for ipair = 1:number_of_pairs
        dev_min(ipair) = min(dev(:,ipair));
        devControl_min(ipair) = min(mean(devControl(:,ipair,:),3)); % or mean across 3D--> mean(devControl,3);
        ratio_strength(ipair) = (dev_min(ipair) - mean(dev(:,ipair))) ./ ...
                                (devControl_min(ipair) - mean(mean(devControl(:,ipair,:),3)));
     end
elseif choice_smooth == 1
    dev_min_smoothed = zeros(1, number_of_pairs);
    devControl_min = zeros(1, number_of_pairs);
    ratio_strength_smooth = zeros(1, number_of_pairs);
    time = (0:150);
    dev_smoothed = zeros(length(time), number_of_pairs);
     for ipair = 1:number_of_pairs
        [dev_smoothed_coeff, S, mu] = polyfit(time, dev(:, ipair)',6);
        dev_smoothed(:, ipair) = polyval(dev_smoothed_coeff, time, S, mu);
        dev_min_smoothed(ipair) = min(dev_smoothed(:, ipair));
        devControl_min(ipair) = min(mean(devControl(:,ipair,:),3)); % or mean across 3D--> mean(devControl,3);
        ratio_strength_smooth(ipair) = (dev_min_smoothed(ipair) - mean(dev_smoothed(:,ipair))) ./ ...
                                (devControl_min(ipair) - mean(mean(devControl(:,ipair,:),3)));
     end    
     ratio_strength = ratio_strength_smooth;
     dev_min = dev_min_smoothed;
     clear dev;
     %dev = dev_smoothed;
end

if strcmp(data_folder, 'pp_batch')
    iweak_pair = find(ratio_strength < 3.50)
    istrong_pair = find(ratio_strength >= 3.50)
elseif strcmp(data_folder, 'pp_poisson')
    iweak_pair = find(ratio_strength < 350)
    istrong_pair = find(ratio_strength >= 350)
end


end