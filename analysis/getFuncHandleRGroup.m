% Look up the handle for a function which works on a set of recordings.
%
% The function handle is taken from a subfunction within this .m file.
% All functions returned have the format
%     X = func(signal, tt, auxlab, ped, vel, lagTime)
% for stimulus-independent or
%     [X, ustm] = func(signal, tt, auxlab, ped, vel, lagTime, stimTimes, ustm)
% for stimulus dependent functions. In the later case, the size of X in
% dimension 2 should match the number of elements in ustm, since this is
% the dimension assumed to denote different stimuli. Alternatively the
% function can return ustm as an empty array if X is not organized by
% stimulus identity.
%
% All inputs should be arrays sized [1, numTimePoints, numRecs], except
% `signal` which should be [numROIs, numTimePoints, numRecs] and both
% `stimTimes` and `ustm`, which will be cell arrays sized [numRecs, 1].
% Here numROIs is the number of non-trivial ROIs present in the recordings,
% numTimePoints is the number of samples over the duration of any single
% recording (the duration in seconds times the sampling frequency), and
% numRecs is the number of recordings to analyse together in the set. This
% input schema means the recordings must all be from the same field of view
% (so they have the same ROIs present) and must have the same duration and
% sampling frequency.
% The output X is sized [numROIs, ...], with the size of higher dimensions
% specific to the function called.
% If there is no function named funcstr within getFuncHandleRGroup.m, but
% there is one in getFuncHandleRec.m, the returned function handle is
% to an anonymous function which applies the single-recording function to
% every recording in the input and concatenates the results along the 3rd
% dimension.
%
% Inputs
% ------
% funcstr : string
%     Name of the function to lookup.
%
% Outputs
% -------
% func : function handle
%     Handle to the requested function.
%
% See also getFuncHandleRec.
%
% ======================================================================= %
%                           LIST OF FUNCTIONS                             %
% ======================================================================= %
% (Ctrl+D to jump to function in this file)
%__________________________________________________________________________
%   # Functions
%   OriSelective
%   OriVar
%   ML_GTdecoder
%   stimDiscrim
%   rwdOnsetData
%   rwdResponsive
%   successRateSMI
%   CorrResponsive
%   CorrSelective
%   tempMatch_taskOri
%   stimDiscrimVR
%   LMI
%   runTime 
%______________________________________________________________________


% ======================================================================= %
%                           MAIN FUNCTION                                 %
% ======================================================================= %

function func = getFuncHandleRGroup(funcstr)

% Look up the function by name from the local scope
% We find the function with the correct name from below and return a handle
% to it!
% But take caution here; str2func always succeeds, even if the function
% doesn't exist.
func = str2func(funcstr);
try
    func(); % This is expected to break due to lack of input
    % If it did not break, raise a warning
    warning('ROCHLAB:functionWithoutInputs', ...
        'Function %s ran without any inputs', funcstr);
catch ME
    if ~strcmp(ME.identifier,'MATLAB:UndefinedFunction') && ...
            ~strcmp(ME.message, ['Undefined function or variable ''' funcstr ''''])
        % If we errored for anything other than a missing function,
        % the function is defined and we have errored due to lack of inputs
        % This is safe to return
        disp('Got function from rec group function handles');
        return;
    end
    % If the function does not exist, try to take it from the single rec
    % functions and stack the outputs together
    try
        innerfunc = getFuncHandleRec(funcstr);
        func = @(varargin) single2group(innerfunc, varargin{:});
        disp('Got function from single rec function handles, and put it in a wrapper');
        warning('ROCHLAB:unnecessarySingleFuncWrap', ...
            'Please consider using calcRecFuncPool with cat instead');
    catch ME2
        if strcmp(ME2.identifier,'MATLAB:UndefinedFunction') && ...
            strcmp(ME2.message, ['Undefined function or variable ''' funcstr '''']);
            rethrow(ME);
        else
            rethrow(ME2);
        end
    end
end

end
function out = single2group(innerfunc, varargin)

% Define dimensions
dROI  = 1;
dTime = 2;
dRec  = 3;

% Check which input args are one-for-every-rec
nArg = length(varargin);
nRecPerArg = cellfun('size', varargin, dRec);
nRec = max(nRecPerArg);

% Loop over every rec
for iRec=1:nRec
    % Assemble the subset of data to put into the single rec function
    args = cell(1,nArg);
    for iArg=1:nArg
        if nRecPerArg(iArg)==1
            args{iArg} = varargin{iArg};
        else
            args{iArg} = varargin{iArg}(:,:,iRec);
        end
    end
    % Apply the inner function to the rec
    X = innerfunc(args{:});
    % Should add some handling for X being a cell array of many outputs...
    % If this is the first rec, we need to initialise the output
    if iRec==1
        if ndims(X)>3
            error('Output has too many dimensions');
        end
        out = zeros(size(X,1),size(X,2),0);
    end
    % Stack results from each rec together
    out = cat(dRec, out, X);
end

end


% ======================================================================= %
%                               FUNCTIONS                                 %
% ======================================================================= %
% Description of possible inputs
%__________________________________________________________________________
% Inputs
% ------
%   signal : array [nROIs, nSamples]
%       Fluorescence activity traces per ROI for each trial and track
%       nSamples = tot number of samples during the recording
%   tt : vector [1, nSamples]
%       Time of samples since the start of the recording (seconds)
%   auxlab -> here, VRlab : structure with fields:
%       Each field is a vector of size [1, nSamples]
%       - trackNum (index of VR track presented by time/frame)
%       - reward (event detection of reward; 0= no reward, 1 = early reward, 2 = default (late) reward)
%       - trialNum (running index of trial number by time/frame)
%       - licks (even detection of licking from lick sensor; 0 = no lick, 1 = lick)
%       - actlab (index of 'action' of animal by time/frame; 0 = stationary; 3 = locomoting)
%   disp : vector [1, nSamples]
%       Cumulative total distance travelled (displacement) since the start 
%       of the recording (cm)
%   vel : vector [1, nSamples]
%       Instantaneous velocity (cm/s)
%   lagTime : vector [1, nSamples]
%       Instantaneous frame by frame lag in VR (seconds)
%__________________________________________________________________________

% ======================================================================= %
%                       FUNCTIONS USING Pre/Post data                     %
% ======================================================================= %

function [X,stimFlag] = OriSelective(signal, tt, auxlab, ~, ~, ~, stimTimes, ustm)
%  Calculates parameters of tuning curve and designates neuron as 
%  Orientation selective based on peak magnitude of the tuning curve vector 
%  (greater than 25th percentile of population) as well as statistical 
%  significance of mean response activity to the prefered vs the orthogonally 
%  oriented grating across trials.     
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 5]
%       Columns:
%           1 - peak angle of tuning curve
%           2 - peak magnitude of the tuning curve
%           3 - prefered stimulus (maximal response)
%           4 - p-value: signed rank test for significance for prefered vs
%               orthogonal oriented grating
%           5 - designation as orientation selective (true ==1)
%
%__________________________________________________________________________

stimFlag = []; 
% Params
oris = [0,45,90,135];
alpha = 0.05;
nROIs = size(signal,1);

% 1) Convert function input to relevant Data 

[dataDf, pref_ori] = prepdataforfunction_RGroup('OriSelective',signal, tt, auxlab, stimTimes, ustm);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%  Data should be matrix [ROIs, Orientations, trials]
%       - dataDf : baseline corrected average DF/F0 for orientations per trial
%  Data should be matrix [ROIs, prefered orientation]
%       - pref_ori : prefered orientations across trials (maximal response)

% initialize output
X = nan(nROIs,5);

% Calculate peak angle and peak magnitude of the tuning curve 
[X(:,1), X(:,2)] = calcOSI(nanmean(dataDf,3), [], oris, true, true); % calculate with circular variance

% get prefered stimulus orientation
X(:,3) = pref_ori;
orth_ori = pref_ori+90;
orth_ori(orth_ori>135)=orth_ori(orth_ori>135)-180;

% get data from prefered and orthogonal orientations
pref_data_index = repmat(oris,nROIs,1)==pref_ori;
orth_data_index = repmat(oris,nROIs,1)==orth_ori;

pref_data = reshape(dataDf(repmat(pref_data_index,1,1,size(dataDf,3))), nROIs, []);
orth_data = reshape(dataDf(repmat(orth_data_index,1,1,size(dataDf,3))), nROIs, []);

% signed rank test for significance for prefered vs orthogonal oriented grating
for iROI = 1:nROIs
    X(iROI,4) = signrank(pref_data(iROI,:), orth_data(iROI,:));
end

% Orientative selective neurons
X(:,5) = (X(:,4) < alpha) & (X(:,2) > prctile(X(:,2),25,1));

end

function [X,stimFlag] = OriVar(signal, tt, auxlab, ~, ~, ~, stimTimes, ustm)
%  Calculates coefficient of variation of the peak magnitude of the tuning 
%  curve across trials for each neuron
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 1]
%       Coefficient of variation across trials (std/mean)
%
%__________________________________________________________________________

stimFlag = []; 
% Params
oris = [0,45,90,135];
nROIs = size(signal,1);

% 1) Convert function input to relevant Data 

dataDf = prepdataforfunction_RGroup('OriSelective',signal, tt, auxlab, stimTimes, ustm);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%  Data should be matrix [ROIs, Orientations, trials]
%       - dataDf : baseline corrected average DF/F0 for orientations per trial

nTrials = size(dataDf,3); % get total number of trials

vector_mag = nan(nROIs, nTrials);
% Calculate peak magnitude of the tuning curve per trial
for iTrial = 1:nTrials
    [~, vector_mag(:,iTrial)] = calcOSI(dataDf(:,:,iTrial), [], oris, true, true); % calculate with circular variance
end

% coefficient of variation across trials (std/mean)
X = nanstd(vector_mag,1,2)./nanmean(vector_mag,2);

end

function [X,stimFlag] = ML_GTdecoder(signal, tt, auxlab, ~, ~, ~, stimTimes, ustm)
%  Bayesian maximum-likelihood decoder 
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, orientation]
%       Decoder accuracy across trials for each orientation
%
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

dataDf = prepdataforfunction_RGroup('OriSelective',signal, tt, auxlab, stimTimes, ustm);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%  Data should be matrix [ROIs, Orientations, trials]
%       - dataDf : baseline corrected average DF/F0 for orientations per trial

[~, X] = run_simple_neuron_ML_GTdecoder(dataDf);

end

function [X,stimFlag] = stimDiscrim(signal, tt, auxlab, ~, ~, ~, stimTimes, ustm)
%  Stimulus discriminability between vertical and angled oriented gratings
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 1]
%       First column: stimulus descriminability (d prime, [d']) between 
%           vertical and angled orientations for each ROI
%
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

[dataDfVert, dataDfAngled] = prepdataforfunction_RGroup('stimDiscrim',signal, tt, auxlab, stimTimes, ustm);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%  Data should be matrix [ROIs, trials]
%       - dataDfVert : average DF/F0 for vertical oriented gratings
%       - dataDfAngled : average DF/F0 for angled oriented gratings

mean1 = nanmean(dataDfVert,2);
mean2 = nanmean(dataDfAngled,2);
var1 = nanvar(dataDfVert,0,2);
var2 = nanvar(dataDfAngled,0,2);

% caclulate d': difference in the means divided by their root mean squared s.d. (variance)
X = (mean1 - mean2) ./ sqrt(0.5 .* (var1 + var2));

end

% ======================================================================= %
%                       FUNCTIONS USING VR data                            %
% ======================================================================= %

function [X,stimFlag] = rwdOnsetData(signal, tt, auxlab, ~, ~, ~, ~, ~)
% Activity surrounding the onset of reward averaged over trials
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nSamples_RewardWindow, nTracks, 2]
%       For each ROI, returns the mean and sem of the average DF/F0 
%       surrounding the onset of reward, for each
%       track or average over all tracks (=>nTracks+1).
%       Dimensions: 
%           1 - ROI
%           2 - samples included in time window around reward
%           3 - VR corridors (tracks)
%           4 - mean [1] or sem [2]
%
%__________________________________________________________________________

stimFlag = [];

% 1) Convert function input to relevant Data 

[ dataDf_RewardWindow, window_samples, window_time, nTracks ] = prepdataforfunction_RGroup('rwdOnsetData',signal, tt, auxlab);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%   dataDf: cell array 
%       A cell can be empty if there was no reward.
%       Otherwise, each cell contains a cell array [1, nTracks], where
%           for each track the DF/F0 of each ROI was extracted for a time 
%           window (2*window_samples) anytime a reward event was detected, 
%           during all trials of the track (nSamples_event)
%           => array [nROIs, 2*window_samples, nSamples_event] 
%   window_samples*2 : number of samples included to extract the signal 
%       surrounding the onset of reward  
%   window_time : time corresponding to window_samples
%   nTracks : number of tracks
%       
% More Parameters
% ===============
f_resample = 40; % re-sampling frequency to equalize datapoints across groups
nSamples_RewardWindow = f_resample*4; % (in frames): get 2 sec before and after Rwd Onset

nROIs = size(signal,1);

% Initialise output
X = nan(nROIs, nSamples_RewardWindow, nTracks, 2);
X_raw = nan(nROIs, window_samples*2, nTracks, 2);

% Compute mean and sem for each track for each ROI
% (signal from onset of reward averaged over trials)
for iTrack = 1:nTracks
    dataDf_track = dataDf_RewardWindow{iTrack}; % size: [nROIs, 2*window_samples, nSamples_event]
    if ~isempty(dataDf_track)
        X_raw(:,:,iTrack,1) = nanmean(dataDf_track,3);
        X_raw(:,:,iTrack,2) = nansem(dataDf_track,3);
    end
end

% Resample to 40Hz
for iTrack = 1:nTracks
    for iData = 1:2
        
        if any(~isnan(X_raw(:,:,iTrack,iData)))
            
            data_resampled = ( resample(X_raw(:,:,iTrack,iData)',window_time,f_resample) )';
            nSamples_resamp = size(data_resampled,2);
            
            % Start filling resampled data
            X(:, 1:min(nSamples_RewardWindow, nSamples_resamp), iTrack, iData) = data_resampled(:, 1:min(nSamples_RewardWindow, nSamples_resamp));
            
            % Check if data is missing
            if nSamples_resamp < nSamples_RewardWindow % pad with last datapoint
                numPad = nSamples_RewardWindow - nSamples_resamp;
                X(:, nSamples_resamp+1:end, iTrack, iData) = repmat(data_resampled(:,end), [1, numPad]);
            end
        end
    end
end

% Replace first/last 5 data points as there can be artifacts from downsampling
X(:,end-5:end,:,:) = X_raw(:,end-5:end,:,:);
X(:,1:5,:,:) = X_raw(:,1:5,:,:);


end

function [X,stimFlag] = rwdResponsive(signal, tt, auxlab, ~, ~, ~, ~, ~)
% Paired t-test for signifiance of increased activity following reward
% onset across trials
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 4]
%       For each ROI, returns the p-value, mean pre and mean post, and  
%       designation as reward responsive [1] or not [0].
%       columns: 
%           1 - p-value from paired t-test (activity in pre vs post reward 
%                  window for each trial)
%           2 - mean pre-reward window (-2 to -1 seconds pre reward onset)
%           3 - mean post-reward window (0 to 1 seconds post reward onset)
%           4 - reward responsive: true == 1
%
%__________________________________________________________________________

stimFlag = [];

% 1) Convert function input to relevant Data 

dataDf_RewardWindow = prepdataforfunction_RGroup('rwdResponsive',signal, tt, auxlab);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%   dataDf_RewardWindow: matrix [ROI, deltaf window, trials]
%       deltaf data for each ROI x window surrounding reward onset (-2 to
%       +2 seconds) x trials
%       
% More Parameters
% ===============
alpha = 0.05; 
nROIs = size(signal,1);

% Initialise output
X = nan(nROIs, 3);

preData = squeeze(nanmean(dataDf_RewardWindow(:,1:size(dataDf_RewardWindow,2)/4,:),2)); % mean pre per trial
postData = squeeze(nanmean(dataDf_RewardWindow(:,size(dataDf_RewardWindow,2)/2:(size(dataDf_RewardWindow,2)/4)*3,:),2)); % mean post per trial

for iROI = 1:nROIs
    [~, X(iROI,1)] = ttest(preData(iROI,:)',postData(iROI,:)');
    X(iROI,2) = nanmean(preData(iROI,:),2); % mean pre per trial
    X(iROI,3) = nanmean(postData(iROI,:),2); % mean post per trial
end

% reward responsive if significant and mean activity post-reward is higher than pre-reward
X(:,4) = (X(:,1) < alpha) & (X(:,3) > X(:,2)); 

end

function [X,stimFlag] = successRateSMI(signal, ~, auxlab, disp, ~, ~, ~, ~)
% Description 
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nTracks]
%       SMI: spatial modulation index. Calculated as percent success rate over trials
%           divided by the shuffled sucess rate for each VR track
%       Note that all rows are just copies as success rate applies to all ROIs.
%
% see also: Shuffle
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

[ success_rate, dataReward, dataLicks, dataDistance, TracksParameters ] = prepdataforfunction_RGroup('successRateSMI', auxlab, disp);

% 2) Analysis 
%
% Format of Data for analysis
% ===========================
%   success_rate : array [1, numTracks]
%           percent success rate over trials for the track
%           NOTE: if there is no data for the track, the entry is NaN
%   dataReward : cell array [1, numTracks]
%       Each cell contains a cell array [1, numTrials] where numTrials is 
%       the number of trials of the track
%   dataLicks : cell array [1, numTracks]
%       Same as dataReward
%   dataDistance : cell array [1, numTracks]
%       Same as dataReward
%   TracksParameters is a structure array [1, numTracks] with fields:
%       - rewardZone (integer)
%       - rewardDefault (integer)
%       These two are parameters to identify the reward zone where 
%       the animal is supposed to lick (see below)
%
% NOTE ABOUT REWARD DATA: 
%           0 = no reward   
%           1 = reward  
%           2 = default reward
%       
% Parameters for analysis
% =======================
numShuffles = 1000;

numTracks = length(TracksParameters);
shuff_rate = nan(numShuffles, numTracks);

for iTrack = 1:numTracks

    if isnan(success_rate(iTrack))
        continue
    end
    
    % Get params of track
    rz = TracksParameters(iTrack).rewardZone;
    default = TracksParameters(iTrack).rewardDefault;
    numTrials = length(dataReward{iTrack});
    
    hit_shuffle = zeros(numShuffles, numTrials);
    
    for iTrial = 1:numTrials
        
        trial_length = length(dataReward{iTrack}{iTrial});
        
        if any(dataReward{iTrack}{iTrial}==1) % hit trial
            Idx_startRwd = find(dataReward{iTrack}{iTrial}==1,1);
            
        elseif any(dataReward{iTrack}{iTrial}==2) % default trial
            Idx_startRwd = find(dataReward{iTrack}{iTrial}==2,1);
            
        else % No reward
            Idx_startRwd = trial_length;
        end
        
        prerwdLicks = dataLicks{iTrack}{iTrial}(1:Idx_startRwd);
        pad = zeros(1,trial_length-Idx_startRwd);
        tmpdata = [pad, prerwdLicks];
        
        shufflick = Shuffle(repmat(tmpdata,numShuffles,1),2);
        
        selected_track_zone = dataDistance{iTrack}{iTrial}>=rz & dataDistance{iTrack}{iTrial}<default ;
        hit_shuffle(:,iTrial) = double( sum(shufflick(:, selected_track_zone), 2)>0 );
    end
    
    shuff_rate(:,iTrack) = sum(hit_shuffle,2)/numTrials; % percent success rate, num hits, num misses
end

% Average over numShuffles
shuff_rate = nanmean(shuff_rate,1); % percent success rate, num hits, num misses

% Ratio
R = success_rate./shuff_rate;

% Simple copy for Output of main function
nROIs = size(signal,1);
X = repmat(R, nROIs, 1);

end

function [X,stimFlag] = CorrResponsive(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Determines if ROIs are responsive to target VR corridors (tracks) presented,
% (i.e. Vertical AND/OR Angled corridor walls) in comparison to reward zone (black walls),
% and returns statistics and mean values for each track.
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nTracks+1, 4]
%       Second dimension: data for each track plus, nTracks+1 data for
%           responsiveness across any target track in place 4 of 3rd dimension
%       Third dimension: For each ROI for each track, 
%           1 - test pre- versus post- reward zone
%               (Wilcoxon signed rank test for zero median)
%           2 - average DF/F0 pre-reward zone
%           3 - average DF/F0 post-reward zone
%           4 - corridor responsive ROIs == 1 (significant and larger mean
%               activity in corridor)
%
% see also: signrank
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 
   
[RZdatapre, RZdatapost] = prepdataforfunction_RGroup('CorrResponsive',signal, tt, auxlab, disp, vel, lagTime);

% 2) Analysis 
%
% Format of Data for analysis
% ===========================
%  Data should be two cell arrays [numTracks, 1]
%       - RZdatapre : DF/F0 for pre-reward zone
%       - RZdatapost : DF/F0 for post-reward zone
%       Each cell contains an array [nROIs, numTrials] where numTrials is
%       the number of trials for the track. The values DF/F0 were obtained
%       by averaging, for each trial, over selected bins pre- / post- rwd.
%       
% Parameters for analysis
% =======================
alpha = 0.001; % 
targetTrack1 = 1; % Angled
targetTrack2 = 3; % Vertical

% Number of tracks and ROIs
nTracks = size(RZdatapre,1);
nROIs = size(signal,1);

% Initialise output array
X = nan(nROIs, nTracks+1, 4); 

for iTrack = 1:nTracks
    
    if ~isempty(RZdatapre{iTrack})
        
        for iROI = 1:nROIs
            X(iROI,iTrack,1) = signrank(RZdatapre{iTrack}(iROI,:),RZdatapost{iTrack}(iROI,:));
        end
        
        X(:,iTrack,2) = nanmean(RZdatapre{iTrack},2); % Corridor activity
        X(:,iTrack,3) = nanmean(RZdatapost{iTrack},2); % Reward zone activity
        
        % determine if ROI is corridor responsive for individual target tracks
        X(:,iTrack,4) = (X(:,iTrack,1) < alpha) & (X(:,iTrack,2) > X(:,iTrack,3)); % Corridor Responsive == 1
    end
    
end

% determine if ROI is corridor responive for any target track
X(:,iTrack+1,4) = X(:,targetTrack1,4) | X(:,targetTrack2,4); % Corridor Responsive == 1

end

function [X,stimFlag] = CorrSelective(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% Determines if ROIs are selective for one of two specific VR corridors (tracks) presented,
% (i.e. Vertical OR Angled corridor walls) in comparison to reward zone (black walls),
% and returns statistics and mean values for each corridor.
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, nTracks+1, 4]
%       Second dimension: data for each track plus, nTracks+1 data for
%           selectivity for a single target track in place 4 of 3rd dimension
%       Third dimension: For each ROI for each track, 
%           1 - test pre- versus post- reward zone
%               (Wilcoxon signed rank test for zero median)
%           2 - average DF/F0 pre-reward zone
%           3 - average DF/F0 post-reward zone
%           4 - corridor responsive ROIs == 1 (significant and larger mean
%               activity in corridor compared to black reward zone
%
% see also: signrank
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 
   
[RZdatapre, RZdatapost] = prepdataforfunction_RGroup('CorrResponsive',signal, tt, auxlab, disp, vel, lagTime);

% 2) Analysis 
%
% Format of Data for analysis
% ===========================
%  Data should be two cell arrays [numTracks, 1]
%       - RZdatapre : DF/F0 for pre-reward zone
%       - RZdatapost : DF/F0 for post-reward zone
%       Each cell contains an array [nROIs, numTrials] where numTrials is
%       the number of trials for the track. The values DF/F0 were obtained
%       by averaging, for each trial, over selected bins pre- / post- rwd.
%       
% Parameters for analysis
% =======================
alpha = 0.001; % 
targetTrack1 = 1; % Angled
targetTrack2 = 3; % Vertical

% Number of tracks and ROIs
nTracks = size(RZdatapre,1);
nROIs = size(signal,1);

% Initialise output array
X = nan(nROIs, nTracks+1, 4); 

for iTrack = 1:nTracks
    
    if ~isempty(RZdatapre{iTrack})
        
        for iROI = 1:nROIs
            X(iROI,iTrack,1) = signrank(RZdatapre{iTrack}(iROI,:),RZdatapost{iTrack}(iROI,:));
        end
        
        X(:,iTrack,2) = nanmean(RZdatapre{iTrack},2); % Corridor activity
        X(:,iTrack,3) = nanmean(RZdatapost{iTrack},2); % Reward zone activity
        
        X(:,iTrack,4) = (X(:,iTrack,1) < alpha) & (X(:,iTrack,2) > X(:,iTrack,3)); % Corridor Responsive == 1
    end
    

end

% determine if ROI is corridor selective for one (but not both) target tracks
X(:,iTrack+1,4) = (X(:,targetTrack1,4) & ~X(:,targetTrack2,4)) | (X(:,targetTrack2,4) & ~X(:,targetTrack1,4)); % Corridor Selective == 1

end

function [X,stimFlag] = tempMatch_taskOri(signal, tt, auxlab, disp, vel, lagTime, ~, ~)
% tepmlate matching decoder computes the decoding accuracy of two target orientations
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 2]
%       Note that all rows are just copies as decoder applies to all ROIs.
%       First column: Give the percentage of correctly classified trials, 
%           for two target orientations (i.e. either Vertical or Angled corridor walls), 
%           given the activities.
%       Second column: Chance level.
%
% see also: templatematching
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

[dataDfVert, dataDfAngled] = prepdataforfunction_RGroup('tempMatch_taskOri',signal, tt, auxlab, disp, vel, lagTime);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%  Data should be two cell arrays [numTracks, 1]
%       - dataDfVert : average DF/F0 for vertical corridor trials
%       - dataDfAngled : average DF/F0 for angled corridor trials
%       Each cell contains an array [nROIs, numTrials] where numTrials is
%       the number of vert or angled trials.
     

% Number of ROIs
nROIs = size(signal,1);

% Initialise output
Xraw = nan(1, 2); 

    if ~isempty(dataDfVert)
        
        numTrials_max = max(size(dataDfVert,2), size(dataDfAngled,2));
        
        dataDfPad = nan(nROIs, 2, numTrials_max); 
                
        for iROI = 1:nROIs
            dataTmp = {dataDfVert(iROI,:) , dataDfAngled(iROI,:)};
            dataDfPad(iROI, :, :) = padcat(dataTmp{:});  
        end
        
        Xraw(1,1) = templatematching(dataDfPad); % decoder accuracy
        Xraw(1,2) = 1/2; % chance level of decoder (i.e. determined by how many track inputs)
        
    end
 
X = repmat(Xraw, nROIs, 1);

end

function [X,stimFlag] = stimDiscrimVR(signal, tt, auxlab, ~, ~, ~, ~, ~)
% Stimulus discriminability between vertical and angled VR corridors
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 1]
%       First column: stimulus descriminability (d prime, [d']) between 
%           vertical and angled corridors for each ROI
%
%__________________________________________________________________________

stimFlag = []; 

% 1) Convert function input to relevant Data 

[dataDfVert, dataDfAngled] = prepdataforfunction_RGroup('stimDiscrimVR',signal, tt, auxlab);

% 2) Analysis
%
% Format of Data for analysis
% ===========================
%  Data should be matrix [ROIs, trials]
%       - dataDfVert : average DF/F0 for vertical corridor trials
%       - dataDfAngled : average DF/F0 for angled corridor trials

mean1 = nanmean(dataDfVert,2);
mean2 = nanmean(dataDfAngled,2);
var1 = nanvar(dataDfVert,0,2);
var2 = nanvar(dataDfAngled,0,2);

X = (mean1 - mean2) ./ sqrt(0.5 .* (var1 + var2));

end

% ======================================================================= %
%                       Generic FUNCTIONS                                 %
% ======================================================================= %

function [X,stimFlag] = LMI(signal, tt, auxlab, ~, ~, ~, stimTimes, ustm)
% Computes the locomotion modulation index. Result is an unbiased estimate
% of the LMI over all stimuli.
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, 1]
%       LMI values across all stimulus categories for each ROI
%__________________________________________________________________________

if isstruct(auxlab) % VR cases
    VRlab = auxlab;
    auxlab = VRlab.actlab;
end

if strcmp(ustm, 'VR') % VR cases
    ustm{1} = {};
end

X = compute_lmi(signal, tt, auxlab, stimTimes, ustm, 0);

stimFlag = [];

end

function [X,stimFlag] = runTime(signal, tt, auxlab, ~, ~, ~, stimTimes, ustm)
% Reports the total number of frames assigned to various behavioural states,
% i.e. stationary, positioning, locomoting, and total frames recorded, for 
% each stimulus period defined in stimTimes input. 
% 
%__________________________________________________________________________
% Outputs
% -------
%   X : array [nROIs, stimulus category[n], behavioural state[4]]
%       Dimension 2, Stimulus category, consists of one column per 
%           stimulus presentation in stimTimes
%       Dimension 3, Behavioural state, is size 4:
%           1 - total frames spent stationary
%           2 - total frames spent positioning (i.e. small movements)
%           3 - total frames spent locomoting 
%           4 - total number of frames across all trials and states
%__________________________________________________________________________


if isstruct(auxlab) % VR cases
    auxlab = auxlab.actlab;
end

nStim = numel(stimTimes{1});
% initialize output
X = nan(size(signal,1), nStim, 4, size(signal,3));

% loop over multiple recording sessions/trials 
for iRec = 1:size(signal,3)
    
    for iStim=1:nStim % loop over different stimulus presentations (i.e. 0, 45, 90, 135, oriented gratings)
        stimTimesTmp = stimTimes{iRec};
        stimTimesTmp(end+1) = size(signal,2)*tt(2);
        % Restrict to only during the presentation of this stimulus
        isDuringStim = stimTimesTmp(iStim) <= tt(:,:,iRec) & tt(:,:,iRec) < stimTimesTmp(iStim+1);
        
        % Take average from each action type
        % Consider stationary
        li = isDuringStim & auxlab(:,:,iRec)==0;
        X(:,iStim,1,iRec) = sum(li,2);
        % Consider positioning
        li = isDuringStim & (auxlab(:,:,iRec)==1 | auxlab(:,:,iRec)==2);
        X(:,iStim,2,iRec) = sum(li,2);
        % Consider locomotion
        li = isDuringStim & auxlab(:,:,iRec)==3;
        X(:,iStim,3,iRec) = sum(li,2);
        % Consider total data
        li = isDuringStim;
        X(:,iStim,4,iRec) = sum(li,2);
        
    end
end
% sum total frames according to action for single stim periods across all trials
X =  squeeze(sum(X,4));

% output order of stimuli (column 2)  
stimFlag = ustm{1};

end

