function [MeanData] = extract_meanDf_grating(signal_df, freq_sampling, stim_time, ustm)
% Extracts response to each stimulus. Outputs the best window of  
% mean_deltaf and a local baseline, and the best orientation.
% 
% Inputs:
%   signal            array   [nRois x nFrames x nTrials]
%   action_labels     array   [1 x nFrames x nTrials]
%   stimID            cell    [nTrials x 1]
%   stimtime          cell    [nTrials x 1]
%
% Output:
%   MeanData          structure array [nRois x 1]
%       struct with fields:
%       - best_orient : the best orientation
%       - grating : response to grating
%           struct with fields:
%           - meanDf : average DF during the whole stim period
%           - meanDf_bestwindow : best window of mean DF
%           - meanDf_baseline : local baseline
    
% _________________________________________________________________________
% Params
% _________________________________________________________________________

% NOTE: Response = meanDf - meanDf_baseline
dt_meanDf = 2; % duration of the the best window of mean df for grating (in seconds)
dt_search_postgrey = 0; % duration of grey after the grating allowed to search for the best window (in seconds)
dt_meanDf_baseline = 2; % duration to evaluate the mean df of the baseline prior the best window (in seconds)
dt_search_baseline_bestwin = 1; % duration AFTER the start of the best window allowed to search for the baseline DF
% E.g. if set to 0, and dt_meanDf_baseline=2s then the last window to
% evaluate the baseline is the 2s prior the start of the best window. If
% set to 1, and dt_meanDf_baseline=2s, it means that the last window to
% evaluate the baseline is centered at the start of the best window.

prop_trim = 0.05; % proportion of stim presentations per orientation to remove before evaluate the best orientation (best mean Response)
% NOTE: prop_trim is applied for both min and max values, hence total trimmed = prop_trim x2

% _________________________________________________________________________

% Small function to ensure that the smooth windows are odd (note: this is
% just to simplify compatibility of functions depending on MATLAB version
% -see function 'locfun_smooth_dF')
func_odd_smoothwin = @(w)2*round(w/2)-1;

% Calculate the number of frames used for the mean DF etc...
nframes_meanDf = func_odd_smoothwin(round(freq_sampling*dt_meanDf));
nframes_search_post = round(freq_sampling*dt_search_postgrey);
nframes_baseline = func_odd_smoothwin(round(freq_sampling*dt_meanDf_baseline));
nframes_search_baseline = round(freq_sampling*dt_search_baseline_bestwin);

[nRois , nFrames , nTrials] = size(signal_df);

% Initialise output
MeanData = struct(...
    'best_orient', cell(nRois, 1), ...
    'grating', cell(nRois, 1));

% Get the stim start/stop frames (greys and gratings)
[static_grating_start_frame, ~, ~, drift_grating_end_frame, grey_start_frame, grey_end_frame, ...
    stim_orientation, stim_direction] = locfun_get_stim_info(ustm, stim_time, freq_sampling);
% Note: all are of size [nStim x nTrials]
nStim = size(stim_orientation, 1);

% A little fix in case imaging stopped too soon
grey_end_frame(grey_end_frame>size(signal_df,2))=size(signal_df,2);
drift_grating_end_frame(drift_grating_end_frame>size(signal_df,2))=size(signal_df,2);

% Organise data per trial
df_grating = cell(nStim, nTrials, nRois);
df_grey = cell(nStim+1, nTrials, nRois); % there is one more grey per trial

for iTrial = 1:nTrials
    
    for iStim = 1:nStim
        
        % Get grating and pre-post grey data
        for iRoi = 1:nRois
            df_grating{iStim, iTrial, iRoi} = signal_df(iRoi, static_grating_start_frame(iStim, iTrial):drift_grating_end_frame(iStim, iTrial), iTrial);
            df_grey{iStim, iTrial, iRoi} = signal_df(iRoi, grey_start_frame(iStim, iTrial):grey_end_frame(iStim, iTrial), iTrial);
        end
      
    end
    
    % Same for the last grey
    for iRoi = 1:nRois
        df_grey{iStim+1, iTrial, iRoi} = signal_df(iRoi, grey_start_frame(iStim+1, iTrial):grey_end_frame(iStim+1, iTrial), iTrial);
    end
    
end


% Loop over ROIs to extract meanDf and meanAction
list_orient = unique(stim_orientation(:,1), 'sorted');
nOrient = length(list_orient);
nDirections = 2;

list_fields = {'meanDf_bestwindow', 'meanDf_baseline', 'meanDf'};
nFields = length(list_fields);
arg = cat(1, list_fields, repmat({[]},1,nFields));
roi_meanData_orient_reset = struct(arg{:});
for iField = 1:nFields
    roi_meanData_orient_reset.(list_fields{iField}) = zeros(nOrient, nDirections, nTrials);
end

num_trim = round(prop_trim*nDirections*nTrials);


for iRoi = 1:nRois
    
    % Extract the best window meanDf during grating, and the
    % corresponding baseline meanDf; also extract the action for each
    % period to time
    [roi_meanData, ~, ~] = locfun_extract_Response(df_grating(:,:,iRoi), df_grey(:,:,iRoi), action_lab_grating, action_lab_grey, ...
        nframes_meanDf, nframes_baseline, nframes_search_post, nframes_search_baseline);
    
 
    % Organise per orientation to find best response Orient and Dir
    roi_meanData_orient = roi_meanData_orient_reset;
    
    for iOrient = 1:nOrient
        is_dir1 = stim_direction == list_orient(iOrient);
        is_dir2 = stim_direction == (list_orient(iOrient)+180);
        
        if size(is_dir2,1)==1 && sum(is_dir2)==0 % for SRP cases with only one long ori
            is_dir2 = stim_direction == list_orient(iOrient);
        elseif sum(is_dir1(:))==numel(stim_direction) % for SRP cases with only one ori
            is_dir1 = false(size(is_dir1));
            is_dir1(1,:) = 1;
            is_dir2 = false(size(is_dir2));
            is_dir2(3,:) = 1;
        end
        % Or could do:
        % is_orient = stim_orientation == list_orient(iOrient);
        % is_dir2 = is_orient & ~is_dir1;
        
        
        
        for iField = 1:nFields
            % Organise per orientation / direction for the whole grating
            % period (static + drift)
            roi_meanData_orient.(list_fields{iField})(iOrient,1,:) = reshape(roi_meanData.(list_fields{iField})(is_dir1), [1, 1, nTrials]);
            roi_meanData_orient.(list_fields{iField})(iOrient,2,:) = reshape(roi_meanData.(list_fields{iField})(is_dir2), [1, 1, nTrials]);

        end
    end
    MeanData(iRoi).grating = roi_meanData_orient;
    
    
    % FIND BEST ORIENTATION USING THE FULL GRATING (STATIC+DRIFT)
    % ** Trim min/max to find best orientation
    R_orient = roi_meanData_orient.meanDf_bestwindow - roi_meanData_orient.meanDf_baseline;
    R_orient_trim = sort(R_orient(:,:), 2, 'descend');
    R_orient_trim = R_orient_trim(:,num_trim+1:end-num_trim);
    % Now find best orientation
    [~, idx_max_orient] = max(mean(R_orient_trim, 2));
    
    MeanData(iRoi).best_orient = list_orient(idx_max_orient);
    

end


end


function [MeanData, idx_start_best_window, idx_end_baseline_allStim, nFrames_stim] = locfun_extract_Response(signal_Df, signal_Df_prepost, rawAction, rawAction_prepost, ...
    nframes_meanDf, nframes_baseline, nframes_search_post, nframes_search_baseline)

search_poststim = nframes_search_post>0;

[nStim, nTrials]= size(rawAction);

% Get the min number of frames during the main stim
nFrames_stim = min(min(cellfun(@(x)length(x), signal_Df)));
nFrames_meanDf_search = nFrames_stim + nframes_search_post;

% PREP PRE/POST STIM
if search_poststim
    
    signal_Df_pre = signal_Df_prepost(1:end-1,:);
    rawAction_pre = rawAction_prepost(1:end-1,:);
    
    signal_Df_post = signal_Df_prepost(2:end,:);
    rawAction_post = rawAction_prepost(2:end,:);
    
else % only care about the signal preceding the DF/F0 of interest
    npre = size(signal_Df_prepost,1);
    if npre > nStim
        signal_Df_pre = signal_Df_prepost(1:end-1,:);
        rawAction_pre = rawAction_prepost(1:end-1,:);
    else
        signal_Df_pre = signal_Df_prepost;
        rawAction_pre = rawAction_prepost;
    end
end



% PREP DATA STIM PERIOD
df_stim_and_post = zeros(nFrames_meanDf_search, nStim, nTrials);
rawAction_stim_and_post = zeros(nFrames_meanDf_search, nStim, nTrials);

for iTrial = 1:nTrials
    for iStim = 1:nStim
        
        % NOTE: take the stim from the end, so make sure get the best part
        if search_poststim
            
            temp_df_stim_and_post = cat(2, signal_Df{iStim, iTrial}(end-nFrames_stim+1:end), signal_Df_post{iStim, iTrial});
            
            temp_rawAction_stim_and_post = cat(2, rawAction{iStim, iTrial}(end-nFrames_stim+1:end), rawAction_post{iStim, iTrial});
            
        else
            
            temp_df_stim_and_post = signal_Df{iStim, iTrial}(end-nFrames_stim+1:end);
            
            temp_rawAction_stim_and_post = rawAction{iStim, iTrial}(end-nFrames_stim+1:end);
        end
        
        df_stim_and_post(:,iStim,iTrial) = temp_df_stim_and_post(1:nFrames_meanDf_search);
        
        rawAction_stim_and_post(:,iStim,iTrial) = temp_rawAction_stim_and_post(1:nFrames_meanDf_search);
        
    end
end


% FIND BEST WINDOW
temp_avgDf = mean(df_stim_and_post(:,:),2); % average all signals
temp_avgDf_movmean = locfun_smooth_dF(temp_avgDf, nframes_meanDf);
[~, idx_start_best_window] = max(temp_avgDf_movmean); % idx_best_window: outputs frame count where the best window starts


% GET RESPONSES and corresponding "BASELINE"(prior best window)
MeanData.meanDf_bestwindow = zeros(nStim, nTrials);
MeanData.meanDf_baseline = zeros(nStim, nTrials);
MeanData.meanDf = zeros(nStim, nTrials);

MeanData.meanAction_bestwindow = struct('still', cell(nStim, nTrials), 'loco', cell(nStim, nTrials));
MeanData.meanAction_baseline = MeanData.meanAction_bestwindow;
MeanData.meanAction = MeanData.meanAction_bestwindow;

idx_end_best_window = idx_start_best_window + nframes_meanDf - 1;
idx_end_baseline_allStim = zeros(nStim, nTrials);


% End frame to search for baseline DF
idx_end_search_baseline = idx_start_best_window + nframes_search_baseline;


for iTrial = 1:nTrials
    for iStim = 1:nStim
        
        % RESPONSE - whole stim signal
        % Save the meanDf over the whole stim (+post if included)
        MeanData.meanDf(iStim, iTrial) = mean(df_stim_and_post(:, iStim, iTrial));
        % Save the action over the whole stim (+post if included)
        MeanData.meanAction(iStim, iTrial) = locfun_get_action(rawAction_stim_and_post(:, iStim, iTrial));
        
        % RESPONSE - Simply take the mean at the best window
        MeanData.meanDf_bestwindow(iStim, iTrial) = mean(df_stim_and_post(idx_start_best_window:idx_end_best_window, iStim, iTrial));
        % Action during stim
        action_labels = rawAction_stim_and_post(idx_start_best_window:idx_end_best_window, iStim, iTrial);
        MeanData.meanAction_bestwindow(iStim, iTrial) = locfun_get_action(action_labels);
        
        
        % BASELINE
        % Prep signal prior the best stim window
        prestim = signal_Df_pre{iStim, iTrial}(end-nframes_baseline+1:end);
        temp_df_before_best_window = cat(1, prestim', df_stim_and_post(1:idx_end_search_baseline,iStim,iTrial));
        % Localise the smallest meanDf window
        temp_meanDf_before_best = locfun_smooth_dF(temp_df_before_best_window, nframes_baseline);
        [MeanData.meanDf_baseline(iStim, iTrial), idx_start_baseline] = min(temp_meanDf_before_best);
        % Action during baseline
        prestim_rawAction = rawAction_pre{iStim, iTrial}(end-nframes_baseline+1:end);
        temp_rawAction_before_best_window = cat(1, prestim_rawAction', rawAction_stim_and_post(1:idx_end_search_baseline,iStim,iTrial));
        action_labels = temp_rawAction_before_best_window(idx_start_baseline:(idx_start_baseline+nframes_baseline-1));
        MeanData.meanAction_baseline(iStim, iTrial) = locfun_get_action(action_labels);
        
        idx_end_baseline_allStim(iStim, iTrial) = idx_start_baseline + nframes_baseline - 1; %output the idx for end baseline regardless of preferred orientation
        
    end
end




end




function [static_grating_start_frame, drift_grating_start_frame, static_grating_end_frame, drift_grating_end_frame, ...
    grey_start_frame, grey_end_frame, ...
    stim_orientation, stim_direction] = locfun_get_stim_info(ustm, stim_time, freq_sampling)

nTrials = length(ustm);

ustm = [ustm{:}]; % cell [nStimuli (black, grey, grating) x nTrials]
stim_time = [stim_time{:}];

all_stim_start_frame = round(stim_time*freq_sampling) + 1; % Don't forget that time 0 is first frame!

% Find the first frame of grating stimuli
is_static = cellfun(@(c)c(1)=='S', ustm);
is_drift = cellfun(@(c)c(1)=='F', ustm);

static_grating_start_frame = reshape(all_stim_start_frame(is_static), [], nTrials);
drift_grating_start_frame = reshape(all_stim_start_frame(is_drift), [], nTrials);

% Find the first frame of grey stimuli
if length(ustm)==364 % bad fix for long familiar trials for active passive experiment
    is_grey = contains(ustm,'b');
    is_grey(1)=0;
    is_grey(end)=0;
else
    try
        is_grey = contains(ustm,'g');
    catch
        is_grey = cellfun(@(c)any(c=='g'), ustm);
    end
end
grey_start_frame = reshape(all_stim_start_frame(is_grey), [], nTrials);

% Find the end frames (ASSUMING GREY, STATIC, DRIFT, GREY)
static_grating_end_frame = drift_grating_start_frame-1;
if size(grey_start_frame,1)<size(drift_grating_start_frame,1) % bad fix for VR active cases with no grey in between protocol
    drift_grating_end_frame = drift_grating_start_frame+(freq_sampling*2); % 2 second drift
else
    drift_grating_end_frame = grey_start_frame(2:end,:)-1;
end
idx_all_greys = reshape(find(is_grey), [], nTrials);
idx_after_last_grey = idx_all_greys(end,:)+1;
after_last_grey_start_frame = reshape(all_stim_start_frame(idx_after_last_grey), [], nTrials);
if size(grey_start_frame,1)<size(drift_grating_start_frame,1) % bad fix for VR active cases with no grey in between protocol
    grey_start_frame = [static_grating_start_frame; grey_start_frame(end,:)-(freq_sampling)*3]; % make grey data equivalent to static 'i.e. pre-drift period')
    grey_end_frame = [static_grating_end_frame; grey_start_frame(end,:)+(freq_sampling)*2]; % make grey data 1 second at end)
elseif isempty(static_grating_start_frame) % fix for SRP case with weird protocol (no Static)
    static_grating_start_frame = static_grating_end_frame-1;
    grey_end_frame = [static_grating_start_frame-1 ; after_last_grey_start_frame-1];
else
    grey_end_frame = [static_grating_start_frame-1 ; after_last_grey_start_frame-1];
end

% Extract angle ID- orientation and direction
static_grating_stim_id = ustm(is_static);
drift_grating_stim_id = ustm(is_drift);

nStimGratings = length(drift_grating_stim_id); %just for size. static and drift same size
stim_orientation = zeros(nStimGratings, 1);
stim_direction = zeros(nStimGratings, 1);
for iGraStim=1:nStimGratings %just for size. static and drift same size\
    if isempty(static_grating_stim_id)
        stim_orientation(iGraStim) = str2double(drift_grating_stim_id{iGraStim}(3:5)); % fix for SRP data missing static
    else
        stim_orientation(iGraStim) = str2double(static_grating_stim_id{iGraStim}(3:5));
    end
    stim_direction(iGraStim) = str2double(drift_grating_stim_id{iGraStim}(3:5));
end
stim_orientation = reshape(stim_orientation, [], nTrials);
if isempty(static_grating_stim_id)
    stim_orientation(stim_orientation>135) = stim_orientation(stim_orientation>135)-180;
end
stim_direction = reshape(stim_direction, [], nTrials);

end


function [dFsmooth] = locfun_smooth_dF(dF, smooth_window)

% NOTE: we know smooth_window is odd! (see start of the script)

try
    dFsmooth = movmean(dF,smooth_window,'Endpoints','discard');
    
catch
    dFsmooth = smooth(dF,smooth_window, 'moving');
    % Need to discard the endpoints...
    left_right_rmv = (smooth_window-1)/2;
    dFsmooth = dFsmooth(left_right_rmv+1:end-left_right_rmv);
end

end

