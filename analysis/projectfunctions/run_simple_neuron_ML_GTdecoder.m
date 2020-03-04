function [DecoderPerf_dir, DecoderPerf_ori] = run_simple_neuron_ML_GTdecoder(GT_R)
% Maximum Likelihood Decoder at the single-neuron level. Orientation decoding is 
% inferred from the direction. Performance is evaluated using leave-one-out procedure:
% the decoder is trained with N-1 trials, and tested on the remaining trial.
%
% Inputs
%   GT_R : array of responses [nRois x nAngles x nTrials]
%
% Outputs
%   DecoderPerf_dir : array [1 x nAngles]
%   DecoderPerf_ori : array [1 x nOrient]


[nRois, nAngles, nTrials] = size(GT_R);
nOrient = nAngles/2;

DecoderPerf_dir = zeros(nRois, nAngles);
DecoderPerf_ori = zeros(nRois, nOrient);

% Preperare the indices of the correct angles
correct_Angle_index = 1:nAngles; % correct direction
correct_Angle180_index = [(nOrient+1):nAngles , 1:nOrient]; % for the orienation, both directions are correct

% List of trial indicices for the leave-one-out procedure
list_idxtrials = 1:nTrials;

%  We leave one out, hence the number of tests is = nTrials
for iTrial = 1:nTrials

    % Data to train the decoder
    X_angle_deco = GT_R(:,:,list_idxtrials(~ (list_idxtrials == iTrial) )); % [nRois x nAngles x (nTrials-1)]
    % Test data
    X_angle_test = GT_R(:,:,iTrial); % [nRois x nAngles x 1]
    
    % Repmat test data
    X_angle_test = repmat(X_angle_test, [1,1,nAngles]); % [nRois x nAngles(test) x nAngles]
        
    % Compute parameters of decoder - approx Gaussian
    decoder_mean = mean(X_angle_deco, 3); % [nRois x nAngles x 1]
    decoder_std = std(X_angle_deco, 0, 3); % [nRois x nAngles x 1]
    
    % Permute and reshape for the next step...
    decoder_mean = permute(decoder_mean, [1,3,2]); % [nRois x 1 x nAngles(deco)]
    decoder_std = permute(decoder_std, [1,3,2]); % [nRois x 1 x nAngles(deco)]
        
    % Likelihood of test data
    % > Transform to zscores for each angle before sending to 'normpdf'
    Z_X_angle_test = (X_angle_test - repmat(decoder_mean, [1,nAngles,1])) ./...
        repmat(decoder_std, [1,nAngles,1]);
    % [nRois x nAngles(test) x nAngles(deco)]
    rois_angle_likelihood = normpdf(Z_X_angle_test);
    % Note: this method works better than below:
%     rois_angle_likelihood = normpdf(X_angle_test, repmat(decoder_mean, [1,nAngles,1]), repmat(decoder_std, [1,nAngles,1]));
    
    % Take log to sum over ROIs instead of product of very small values
    angle_log_likelihood = log(rois_angle_likelihood); % [nRois x nAngles(test) x nAngles(deco)]

    % Get MAX argument (i.e. angle that leads to highest likelihood)
    [~, angle_max_index] = max(angle_log_likelihood, [], 3); % [nRois x nAngles(test) x 1]
 
    % >>> Performance per roi
    is_correct_angle = angle_max_index == repmat(correct_Angle_index,[nRois,1]); % [nRois x nAngles(test)]
    is_correct_angle180 = angle_max_index == repmat(correct_Angle180_index,[nRois,1]); % [nRois x nAngles(test)]
    is_correct_orient = is_correct_angle | is_correct_angle180; % for the orienation, both directions are correct
    
    % Save the new performance 
    DecoderPerf_dir = DecoderPerf_dir + double(is_correct_angle);
    % The two directions contribute to the performance of orientation decoding
    DecoderPerf_ori = DecoderPerf_ori + double(is_correct_orient(:,1:nOrient)) + double(is_correct_orient(:,nOrient+1:end));

end

% Finally, compute the average performance
DecoderPerf_ori = DecoderPerf_ori / (2*nTrials); % [nRois x nOrient]
DecoderPerf_dir = DecoderPerf_dir / nTrials; % [nRois x nAngles]

end


