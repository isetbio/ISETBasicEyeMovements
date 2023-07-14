function [responses] = computeConeResponseforSVM(stimParams, sceneFlag, varargin)
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
%    'nTrialsNum'                   - Integer. Specify the number of trials
%                                     created for each time point
%    'center'                       - [x y] of the center for the grating. 
%                                     default [0 0]
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
    p.addParameter('nTrialsNum',10, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('center', [0 0], @(x) (isempty(x) | isscalar(x)));
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    
    % Retrieve the response noiseFlag labels and validate them.
    viewingDistance = p.Results.viewingDistance;
    minimumPixelsPerHalfPeriod = p.Results.minimumPixelsPerHalfPeriod;
    integrationTime = p.Results.integrationTime;
    nTrialsNum = p.Results.nTrialsNum;
    center = p.Results.center;

    % Generate Scene based on stimulus params
    presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', viewingDistance); 
    scene = generateGaborSceneCustomize(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay, ...
        'minimumPixelsPerHalfPeriod', minimumPixelsPerHalfPeriod);
    % Generate human optics
    theOI = oiCreate('wvf human');
    OI = oiCompute(theOI, scene);
    % Create Cone Mosaic
    theMosaic = cMosaic('sizeDegs', [1, 1] * stimParams.sizeDegs, ...
                'integrationTime', integrationTime / 1000);
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
    noisyExcitation = zeros(nTrialsNum,nTimebin*nFrames,3822);
    if strcmp(sceneFlag,'test')
        stimParams.center = center;
        scene_drift = generateGaborSceneCustomize(...
            'stimParams', stimParams,...
            'presentationDisplay', presentationDisplay, ...
            'minimumPixelsPerHalfPeriod', 10);
        OI_drift = oiCompute(theOI, scene_drift);
        stimParams.contrast = 0;
        scene_null = generateGaborSceneCustomize(...
            'stimParams', stimParams,...
            'presentationDisplay', presentationDisplay, ...
            'minimumPixelsPerHalfPeriod', 10);
        OI_null = oiCompute(theOI, scene_null);
        for f = 1:nFrames
            if f == 4
                [noiseFreeExcitation_temp,noisyExcitation_temp,~,~,~] = theMosaic.compute(OI, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
            elseif f == 5
                [noiseFreeExcitation_temp,noisyExcitation_temp,~,~,~] = theMosaic.compute(OI_drift, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
            else
                [noiseFreeExcitation_temp,noisyExcitation_temp,~,~,~] = theMosaic.compute(OI_null, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
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
        [noiseFreeExcitation_temp,noisyExcitation_temp,~,~,~] = theMosaic.compute(OI_null, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
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
    noiseFreePhotocurr = computePhotocurrent(noiseFreeExcitation, timeAxis8frames, 'none');
    noisyPhotocurr = computePhotocurrent(noisyExcitation, timeAxis8frames, 'random');
    responses.noiseFreeExcitation = noiseFreeExcitation;
    responses.noisyExcitation = noisyExcitation;
    responses.noiseFreePhotocurr = noiseFreePhotocurr;
    responses.noisyPhotocurr = noisyPhotocurr;
end
