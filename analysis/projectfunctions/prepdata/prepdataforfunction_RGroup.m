function [varargout] = prepdataforfunction_RGroup(nameOfFunction, varargin)

List_func_do_VRdata2tracks = {...
    'CorrResponsive','tempMatch_taskOri'};


if ismember(nameOfFunction, List_func_do_VRdata2tracks)
    % *********************************************************************
    % -- GROUP OF FUNCTIONS CALLING VRdata2tracks
    % *********************************************************************
    
    signal_df = varargin{1};
    time_samples = varargin{2};
    VRlab = varargin{3};
    disp = varargin{4};
    vel = varargin{5};
    lagTime = varargin{6};
    % DEFAULT PARAMS
    dataTrialsOnly = []; % if true, only returns trials with data (remove nan trials of other tracks or categories)
    useHitsOnly = false; % if true, only successful trials are used (licked in reward zone)
    useMissOnly = false;
    
    
    % Organize by track/distance etc
    [~, dFTrials] = ...
        VRdata2tracks(signal_df, time_samples, VRlab, disp, vel, lagTime, [], [], [], dataTrialsOnly, useHitsOnly, useMissOnly);
    % Note:
    % dFTrials:  cell array [numTracks, 1]
    %   each cell is: array [nROIs, numBins, numTrials]
    %       Binned fluorescence dF/F signal for each trial of the track
    
    % Check empty tracks
    nonemptyTracks = ~ cellfun(@isempty, dFTrials);
    
    % Get number of tracks
    nTracks = length(dFTrials);
    
    % Set output depending on function requested
    switch nameOfFunction
        
        case 'CorrResponsive'
            % Note: RZ starts 16 bins (40cm) from end; 120cm length track RZ starts at bin 32
            % To have same amt of data (10 bins)
            dataRangepre = [30,20]; % from end of trial; Gratings: 10 bins total (25cm), from end [30,20] = -35 to -10 cm from reward zone onset.
            dataRangepost = [16,6]; % from end of trial; Reward Zone: 10 bins total (25cm), from end [16,6] = 0 to +25 cm from reward zone onset.
            
            RZdatapre = cell(nTracks, 1);
            RZdatapost = cell(nTracks, 1);
            
            RZdatapre(nonemptyTracks) = cellfun(@squeeze, (cellfun(@(x) nanmean(x(:,end-dataRangepre(1):end-dataRangepre(2),:),2), dFTrials(nonemptyTracks), 'UniformOutput', false)),'UniformOutput', false);
            RZdatapost(nonemptyTracks) = cellfun(@squeeze, (cellfun(@(x) nanmean(x(:,end-dataRangepost(1):end-dataRangepost(2),:),2), dFTrials(nonemptyTracks), 'UniformOutput', false)),'UniformOutput', false);
            
            varargout{1} = RZdatapre;
            varargout{2} = RZdatapost;
            
        case 'tempMatch_taskOri'
            
            dataDf = cell(nTracks, 1);
            dataDf(nonemptyTracks) = cellfun(@squeeze, (cellfun(@(x) nanmean(x,2), dFTrials(nonemptyTracks), 'UniformOutput', false)),'UniformOutput', false);
            
            % Separate the sucessful (Hit) and miss trials
            dataDfVert = dataDf{3};
            dataDfAngled = dataDf{1};
            
            varargout{1} = dataDfVert;
            varargout{2} = dataDfAngled;
            
    end
    
else
    
    switch nameOfFunction
        
        % *********************************************************************
        % -- VARIOUS FUNCTIONS USING VR data
        % *********************************************************************
        case 'rwdOnsetData'
            
            signal_df = varargin{1};
            time_samples = varargin{2};
            VRlab = varargin{3}; % note: all fields are [1, numSamples]
            
            % DEFAULT PARAMS of tracks
            [numTracks] = get_tracks_info();
            
            avgHz = 1/(mean(diff(time_samples)));
            window_samples = floor(avgHz*2); % (in frames): avgHz = 1000ms window; get 2sec before and after Rwd Onset
            window_time = time_samples(1:window_samples*2);
            nROIs = size(signal_df,1);
            
            numSamples = length(VRlab.trackNum);
            list_indices = 1:numSamples;
            
            % NOTE: skip last trial in case it is incomplete
            is_possible_Index = ( list_indices >= (window_samples*2 +1) ) & ( (list_indices+window_samples) <= numSamples );
            
            % Find reward events
            is_reward_event = VRlab.reward>0;
            
            % Make sure that an event was not too soon or too late
            is_reward_event = is_reward_event & is_possible_Index ;
            list_indices_selected = list_indices(is_reward_event);
            
            numSamples_events = length(list_indices_selected);
            
            eventAvg = nan(nROIs, window_samples*2, numSamples_events);
            for iEvent = 1:numSamples_events
                selected_samples = list_indices_selected(iEvent)-(window_samples-1) : list_indices_selected(iEvent)+window_samples ;
                % Take mean data across trials from window after events (onset | dataWindow)
                eventAvg(:,:,iEvent) = signal_df(:,selected_samples);
            end
            
            dataDf_RewardWindow = cell(1,numTracks);
            for iTrack = 1:numTracks
                is_event_of_track = VRlab.trackNum(list_indices_selected) == iTrack;
                dataDf_RewardWindow{iTrack} = eventAvg(:,:,is_event_of_track);
            end
            
            varargout{1} = dataDf_RewardWindow;
            varargout{2} = window_samples;
            varargout{3} = window_time;
            varargout{4} = numTracks;
            
        case 'rwdResponsive'
            
            signal_df = varargin{1};
            time_samples = varargin{2};
            VRlab = varargin{3}; % note: all fields are [1, numSamples]
            
            % PARAMS
            avgHz = 1/(mean(diff(time_samples)));
            window_samples = floor(avgHz*2); % (in frames): avgHz = 1000ms window; get 2sec before and after Rwd Onset
            nROIs = size(signal_df,1);
            
            numSamples = length(VRlab.trackNum);
            list_indices = 1:numSamples;
            
            % NOTE: skip last trial in case it is incomplete
            is_possible_Index = ( list_indices >= (window_samples*2 +1) ) & ( (list_indices+window_samples) <= numSamples );
            
            % Find reward events
            is_reward_event = VRlab.reward>0;
            
            % Make sure that an event was not too soon or too late
            is_reward_event = is_reward_event & is_possible_Index ;
            list_indices_selected = list_indices(is_reward_event);
            
            numSamples_events = length(list_indices_selected);
            
            dataDf_RewardWindow = nan(nROIs, window_samples*2, numSamples_events);
            for iEvent = 1:numSamples_events
                selected_samples = list_indices_selected(iEvent)-(window_samples-1) : list_indices_selected(iEvent)+window_samples ;
                % Take mean data across trials from window after events (onset | dataWindow)
                dataDf_RewardWindow(:,:,iEvent) = signal_df(:,selected_samples);
            end
            
            varargout{1} = dataDf_RewardWindow;
            
        case 'stimDiscrimVR'
            
            signal_df = varargin{1};
            time_samples = varargin{2};
            VRlab = varargin{3}; % note: all fields are [1, numSamples]
            
            % DEFAULT PARAMS
            targetTrack1 = 3; % Vertical corridor
            targetTrack2 = 1; % Angled corridor
            
            avgHz = 1/(mean(diff(time_samples))); % avgHz = 1 second window;
            window_samples = floor(avgHz*2); % (in frames): get 2 seconds after Trial Onset
            nROIs = size(signal_df,1);
            
            % Find trial start for first target corridor
            list_index = VRlab.trackNum == targetTrack1;
            listTrials = unique(VRlab.trialNum(list_index));
            % Initialize output
            dataDftarget1 = nan(nROIs, length(listTrials));
            
            % get data from time window per trial for target 1 corridor
            for iTrial = 1:length(listTrials)
                is_trial_startWindow = find(VRlab.trialNum==listTrials(iTrial),1);
                % Take mean across window after trial onset for each trial
                dataDftarget1(:, iTrial) = nanmean(signal_df(:, is_trial_startWindow : is_trial_startWindow + window_samples -1),2);
            end
            
            % Find trial start for second target corridor
            list_index = VRlab.trackNum == targetTrack2;
            listTrials = unique(VRlab.trialNum(list_index));
            % Initialize output
            dataDftarget2 = nan(nROIs, length(listTrials));
            
            % get data from time window per trial for target 1 corridor
            for iTrial = 1:length(listTrials)
                is_trial_startWindow = find(VRlab.trialNum==listTrials(iTrial),1);
                % Take mean across window after trial onset for each trial
                dataDftarget2(:, iTrial) = nanmean(signal_df(:, is_trial_startWindow : is_trial_startWindow + window_samples -1),2);
            end
            
            varargout{1} = dataDftarget1;
            varargout{2} = dataDftarget2;
            
        case 'successRateSMI'
            
            VRlab = varargin{1}; % note: all fields are [1, numSamples]
            track_distance = varargin{2};
            
            % DEFAULT PARAMS of tracks
            [numTracks, TracksParameters] = get_tracks_info();
            
            success_rate = nan(1, numTracks);
            dataReward = cell(1,numTracks);
            dataLicks = cell(1,numTracks);
            dataDistance = cell(1,numTracks);
            
            % Calculate success rate and get reward & lick data
            for iTrack = 1:numTracks
                track_samples = VRlab.trackNum == iTrack;
                if ~any(track_samples)
                    continue
                end
                
                list_trialNumbers_track = unique(VRlab.trialNum(track_samples));
                numTrials = length(list_trialNumbers_track);
                dataReward{iTrack} = cell(1, numTrials);
                dataLicks{iTrack} = cell(1, numTrials);
                dataDistance{iTrack} = cell(1, numTrials);
                numhitTrials = 0;
                for iTrial = 1:numTrials
                    selected_trial = list_trialNumbers_track(iTrial);
                    trial_track_samples = track_samples & (VRlab.trialNum == selected_trial);
                    if any(VRlab.reward(trial_track_samples)==1)
                        numhitTrials = numhitTrials+1;
                    end
                    
                    dataReward{iTrack}{iTrial} = VRlab.reward(trial_track_samples);
                    dataLicks{iTrack}{iTrial} = VRlab.licks(trial_track_samples);
                    dataDistance{iTrack}{iTrial} = track_distance(trial_track_samples);
                end
                
                success_rate(iTrack) = numhitTrials/numTrials; % percent success rate, num hits, num misses
            end
            
                       
            % Output
            varargout{1} = success_rate;
            varargout{2} = dataReward;
            varargout{3} = dataLicks;
            varargout{4} = dataDistance;
            varargout{5} = TracksParameters;
            
            
        % *****************************************************************
        % -- VARIOUS FUNCTIONS USING Pre/Post data
        % *****************************************************************
        case 'OriSelective'
            
            signal_df = varargin{1};
            time_samples = varargin{2};
            auxlab = varargin{3};
            stimTimes = varargin{4};
            ustm = varargin{5};
            
            if isstruct(auxlab) % VR cases
                auxlab = auxlab.actlab;
            end
            freq_sampling = round(1/(mean(diff(time_samples))));
            MeanData = extract_meanDf_grating(signal_df, freq_sampling, stimTimes, ustm);

            % DEFAULT PARAMS
            nROIs = size(signal_df,1);
            nTrials = length(stimTimes)*2;
            oris = [0,45,90,135];
            
            dataDf = nan(nROIs,length(oris),nTrials);
            pref_ori = nan(nROIs,1);
            
            for iROI = 1:nROIs
                for iOri = 1:length(oris)
                    % reshape for one trial per orientation presentation [ROIs x trials]
                    dataDf(iROI,iOri,:) = reshape(...
                        MeanData(iROI).grating.meanDf_bestwindow(iOri,:,:) - ...
                        MeanData(iROI).grating.meanDf_baseline(iOri,:,:), 1, []);
                    pref_ori(iROI,1) = MeanData(iROI).best_orient;
                end
            end
            
            varargout{1} = dataDf;
            varargout{2} = pref_ori;
        
        case 'ML_GTdecoder'
            
            signal_df = varargin{1};
            time_samples = varargin{2};
            auxlab = varargin{3};
            stimTimes = varargin{4};
            ustm = varargin{5};
            
            freq_sampling = round(1/(mean(diff(time_samples))));
            MeanData = extract_meanDf_grating(signal_df, freq_sampling, stimTimes, ustm);
            
            % DEFAULT PARAMS
            nROIs = size(signal_df,1);
            nTrials = length(stimTimes)*2;
            oris = [0,45,90,135];
            
            dataDf = nan(nROIs,length(oris),nTrials);
            
            for iROI = 1:nROIs
                for iOri = 1:length(oris)
                    % reshape for one trial per orientation presentation [ROIs x trials]
                    dataDf(iROI,iOri,:) = reshape(MeanData(iROI).grating.meanDf(iOri,:,:), 1, []);
                end
            end
            
            varargout{1} = dataDf;
            
        case 'stimDiscrim'
            
            signal_df = varargin{1};
            time_samples = varargin{2};
            auxlab = varargin{3};
            stimTimes = varargin{4};
            ustm = varargin{5};
            
            freq_sampling = round(1/(mean(diff(time_samples))));
            MeanData = extract_meanDf_grating(signal_df, freq_sampling, stimTimes, ustm);
            
            % DEFAULT PARAMS
            nROIs = size(signal_df,1);
            nTrials = length(stimTimes)*2;
            % oris: [0,45,90,135]
            targetOri1 = 3; % Vertical corridor [90]
            targetOri2 = 2; % Angled corridor [45]
            
            dataDftarget1 = nan(nROIs,nTrials);
            dataDftarget2 = nan(nROIs,nTrials);
            
            for iROI = 1:nROIs
                % reshape for one trial per orientation presentation [ROIs x trials]
                dataDftarget1(iROI,:) = reshape(MeanData(iROI).grating.meanDf(targetOri1,:,:), 1, []);
                dataDftarget2(iROI,:) = reshape(MeanData(iROI).grating.meanDf(targetOri2,:,:), 1, []);
            end
            
            
            varargout{1} = dataDftarget1;
            varargout{2} = dataDftarget2;
    end
    
end



end

%__________________________________________________________________________
% ___ FUNCTIONS ___________________________________________________________


function [numTracks, TracksParameters] = get_tracks_info()


numTracks = 8;
% Tracks parameters
TracksParameters = struct('rewardZone', cell(1,numTracks), 'rewardDefault', cell(1,numTracks));
[TracksParameters([1,2]).rewardZone] = deal(80);
[TracksParameters([1,2]).rewardDefault] = deal(100);
[TracksParameters([3,5]).rewardZone] = deal(120);
[TracksParameters([3,5]).rewardDefault] = deal(140);
[TracksParameters([7,8]).rewardZone] = deal(200);
[TracksParameters([7,8]).rewardDefault] = deal(220);

end
