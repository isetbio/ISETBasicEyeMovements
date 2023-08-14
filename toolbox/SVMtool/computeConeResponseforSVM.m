function [responses] = computeConeResponseforSVM(stimParams, theMosaic, sceneFlag, varargin)

% Compute function for computation of photopigment excitation and
% photocurrent for the corresponding scene
%
% Syntax:
%   responses = computeConeResponsesforSVM(...
%    stimParams, sceneFlag, varargin);
% Inputs:
%    stimParams                     - a struct containing properties of the
%                                     the grating.
%    sceneFlag                      - a string specifiying the type of
%                                     scene that desired to create.
%                                     valid string are:
%                                       - 'null': all 8 scenes are background without any gratings
%                                       - 'test': follow the pattern of 3-2-3 
%                                          (3 are backgrounds and 2 are gratings)
% Optional key/value input arguments:
%    'viewingDistance'              - Double. default = 0.50     
%    'minimumPixelsPerHalfPeriod'   - Integer. default = 10
%    'integrationTime'              - Integer. integration time for cMosaic
%                                     object. default = 15 
%    'nTrials'                   - Integer. Specify the number of trials
%                                     created for each time point
%    'center'                       - [x y] of the center for the grating. 
%                                     default [0 0]
%    'visualize'                    - String. flag to indicate whether to
%                                     visualize the response. Valid String are:
%                                       - 'excitation'
%                                       - 'photocurrent'
%                                       - 'all'
%    'contrast'                     - Double between 0-100. Used for changing the
%                                     contrast of the test stimulus.
%    'responseFlag'                 - Specify expected response type to
%                                     save time. valid string is
%                                     'excitation'. If not specified,
%                                     calculate all responses.
% 
% Outputs:
%    responses  - A struct that contains cone excitation and photocurrent and timeAxis. 
%                 When access specific responses, refers to:
%                 responses.noiseFreeExcitation
%                 responses.noisyExcitation
%                 responses.noiseFreePhotocurr
%                 responses.noisyPhotocurr
%                 responses.timeAxis
%
% History:
%    13/07/23  Qiyuan Feng

    % Parse the input arguments
    p = inputParser;
    p.addParameter('viewingDistance', 0.50, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('minimumPixelsPerHalfPeriod', 10, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('integrationTime',15, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('nTrials',10, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('center', [0 0]);
    p.addParameter('contrast', '', @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('visualize', '', @(x) (isempty(x) | ischar(x)));
    p.addParameter('responseFlag', '', @(x) (isempty(x) | ischar(x)));
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    
    % Retrieve the response noiseFlag labels and validate them.
    viewingDistance = p.Results.viewingDistance;
    minimumPixelsPerHalfPeriod = p.Results.minimumPixelsPerHalfPeriod;
    integrationTime = p.Results.integrationTime;
    nTrials = p.Results.nTrials;
    center = p.Results.center;
    contrast = p.Results.contrast;
    visualize = p.Results.visualize;
    responseFlag = p.Results.responseFlag;

    % Generate Scene based on stimulus params
    presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', viewingDistance);
    if ~isempty(contrast)
        stimParams.contrast = contrast/100;
    end
    % disp(stimParams.contrast)
    scene = generateGaborSceneCustomize(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay, ...
        'minimumPixelsPerHalfPeriod', minimumPixelsPerHalfPeriod);
    % visualizeScene(scene, 'displayRadianceMaps', false);
    % Generate human optics
    theOI = oiCreate('wvf human');
    OI = oiCompute(theOI, scene);
    % Create Cone Mosaic
%     theMosaic = cMosaic('sizeDegs', [1, 1] * stimParams.sizeDegs, ...
%                 'integrationTime', integrationTime / 1000);
    % generate the experiment sequence of 8 frames: 3 - 2 - 3
    nTimebin = integrationTime;
    nFrames = 8;
    ntimeStep = (integrationTime / 1000)/ nTimebin;
    timeAxis8frames = [];
    for t = 0+ntimeStep:ntimeStep:(integrationTime/1000)*nFrames
        timeAxis8frames(end+1) = t;
    end
    responses.timeAxis = timeAxis8frames;
    % cone excitation
    noiseFreeExcitation = zeros(1,nTimebin*nFrames,3822);
    noisyExcitation = zeros(nTrials,nTimebin*nFrames,3822);
    
    if strcmp(sceneFlag,'test')
        stimParams.center = center;
        scene_drift = generateGaborSceneCustomize(...
            'stimParams', stimParams,...
            'presentationDisplay', presentationDisplay, ...
            'minimumPixelsPerHalfPeriod', 10);
        % visualizeScene(scene_drift, 'displayRadianceMaps', false);
        OI_drift = oiCompute(theOI, scene_drift);
        stimParams.contrast = 0;
        scene_null = generateGaborSceneCustomize(...
            'stimParams', stimParams,...
            'presentationDisplay', presentationDisplay, ...
            'minimumPixelsPerHalfPeriod', 10);
        % visualizeScene(scene_null, 'displayRadianceMaps', false);
        OI_null = oiCompute(theOI, scene_null);
        for f = 1:nFrames
            if f == 4
                [noiseFreeExcitation_temp,noisyExcitation_temp,~,~,~] = theMosaic.compute(OI, 'nTrials', nTrials, 'withFixationalEyeMovements', false); % OI
            elseif f == 5
                [noiseFreeExcitation_temp,noisyExcitation_temp,~,~,~] = theMosaic.compute(OI_drift, 'nTrials', nTrials, 'withFixationalEyeMovements', false); % OI_drift
            else
                [noiseFreeExcitation_temp,noisyExcitation_temp,~,~,~] = theMosaic.compute(OI_null, 'nTrials', nTrials, 'withFixationalEyeMovements', false);
            end
            for i = 1:nTimebin
                count = i + nTimebin*(f-1);
                noiseFreeExcitation(:,count,:) = noiseFreeExcitation_temp;
                noisyExcitation(:,count,:) = noisyExcitation_temp;
            end
        end
    elseif strcmp(sceneFlag, 'null')
        % force to set the contrast to 0 just in case that of forgetting to
        % change the contrast when creating null stimulus
        stimParams.contrast = 0;
        scene_null = generateGaborSceneCustomize(...
            'stimParams', stimParams,...
            'presentationDisplay', presentationDisplay, ...
            'minimumPixelsPerHalfPeriod', 10);
        OI_null = oiCompute(theOI, scene_null);
        for f = 1:nFrames
            [noiseFreeExcitation_temp,noisyExcitation_temp,~,~,~] = theMosaic.compute(OI_null, 'nTrials', nTrials, 'withFixationalEyeMovements', false);
            for i = 1:nTimebin
                count = i + nTimebin*(f-1);
                noiseFreeExcitation(:,count,:) = noiseFreeExcitation_temp;
                noisyExcitation(:,count,:) = noisyExcitation_temp;
            end
        end
    else
        error('invalid scene type');
    end
    % photocurrent
    if strcmp(responseFlag, 'excitation')
        noiseFreePhotocurr = NaN;
        noisyPhotocurr = NaN;
    else
        noiseFreePhotocurr = computePhotocurrent(noiseFreeExcitation, timeAxis8frames, 'none');
        noisyPhotocurr = computePhotocurrent(noisyExcitation, timeAxis8frames, 'random');
    end
    
    responses.noiseFreeExcitation = noiseFreeExcitation;
    responses.noisyExcitation = noisyExcitation;
    responses.noiseFreePhotocurr = noiseFreePhotocurr;
    responses.noisyPhotocurr = noisyPhotocurr;
    
    if strcmp(visualize,'excitation')
        visualizeResponseOverTime(responses,visualize);
    elseif strcmp(visualize,'photocurret')
        visualizeResponseOverTime(responses,visualize);
    elseif strcmp(visualize,'all')
        visualizeAllResponses(responses);
    end
end
