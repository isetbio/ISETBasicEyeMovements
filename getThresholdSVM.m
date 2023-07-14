function threshold = getThresholdSVM(varargin)

% Compute function for detection threshold of different grating contrast
%
% Syntax:
%   threshold = getThresholdSVM(varargin);
%
% Inputs:
%    neuralEngineOBJ                - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    neuralResponseParamsStruct     - a struct containing properties of the
%                                     employed neural chain.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%    instancesNum                   - the number of response instances to compute
%
% Optional key/value input arguments:
%    'noiseFlags'                   - Cell array of strings containing labels
%                                     that encode the type of noise to be included
%                                     Valid values are: 
%                                        - 'none' (noise-free responses)
%                                        - 'random' (noisy response instances)
%                                     Default is {'random'}.
%   'rngSeed'                       - Integer.  Set rng seed. Empty (default) means don't touch the
%                                     seed.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams (optics and midget RGC mosaic params) that define the neural 
%               compute pipeline for this computation.  This can be useful
%               for a user interested in knowing what needs to be supplied
%               to this.
%
%             - If called from a parent @neuralResponseEngine), the returned
%               struct is organized as follows:
%                .neuralResponses : dictionary of responses indexed with 
%                                   labels corresponding to the entries of
%                                   the 'noiseFlags'  optional argument
%                .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%                .neuralPipeline  : a struct containing the optics and the mRGC mosaic
%                                   employed in the computation (only returned if 
%                                   the parent @neuralResponseEngine object has 
%                                   an empty neuralPipeline property)
%
%       The computed neural responses can be extracted as:
%           neuralResponses('one of the entries of noiseFlags') 
%       and are arranged in a matrix of:
%           [instancesNum x mCones x tTimeBins] 
%

% History:
%    13/07/23  Qiyuan Feng
    
end
%% Create Stimulus & Response
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 0.50);
% Stimulus 1 
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 10, ... 
    'orientationDegs', 90, ...               
    'phaseDegs', 0, ...                    
    'sizeDegs', 0.5, ...                    
    'sigmaDegs', 0.25/3, ...                 
    'contrast', 100/100,...                  
    'meanLuminanceCdPerM2', 15, ...
    'center', [0 0], ...
    'pixelsAlongWidthDim', [], ...          
    'pixelsAlongHeightDim', [] ...          
    );
scene1 = generateGaborSceneCustomize(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay, ...
    'minimumPixelsPerHalfPeriod', 10);

stimParams.center = [0 0];
scene1_2Frames = generateGaborSceneCustomize(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay, ...
    'minimumPixelsPerHalfPeriod', 10);

% Stimulus 2: NULL
stimParams.contrast = 0;
scene2 = generateGaborSceneCustomize(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay, ...
    'minimumPixelsPerHalfPeriod', 10);
    % visualizeScene(scene1, ...
    %     'displayRadianceMaps', false);.
    % visualizeScene(scene2, ...
    %     'displayRadianceMaps', false);

% Generate human optics
theOI = oiCreate('wvf human');
OI_1 = oiCompute(theOI, scene1);
OI_2 = oiCompute(theOI, scene2);
OI_1_2Frames = oiCompute(theOI, scene1_2Frames);
    % visualizeOpticalImage(OI_1, ...
    %     'displayRadianceMaps', false);
    % visualizeOpticalImage(OI_2, ...
    %     'displayRadianceMaps', false);

% Create Cone Mosaic
integrationTime = 15;
theMosaic = cMosaic('sizeDegs', [1, 1] * stimParams.sizeDegs, ...
            'integrationTime', integrationTime / 1000);

%% 8 Frames Experiment: 3 - 2 - 3
nTrialsNum = 10;
nTimebin = 15;
nFrames = 8;

ntimeStep = (integrationTime / 1000)/ nTimebin;
timeAxis8frames = [];
for t = 0+ntimeStep:ntimeStep:(integrationTime/1000)*nFrames
    timeAxis8frames(end+1) = t;
end
% Cone Excitation
% Null sequence (8 background frames)
[noiseFreeExcitationNull_same, noiseFreeExcitationTest_same] = deal(zeros(1,nTimebin*nFrames,3822));
[noisyExcitationNull_same, noisyExcitationTest_same] = deal(zeros(nTrialsNum,nTimebin*nFrames,3822));
for f = 1:nFrames
    [noiseFreeExcitationNull_temp,noisyExcitationNull_temp,~,~,~] = theMosaic.compute(OI_2, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
    for i = 1:nTimebin
        count = i + nTimebin*(f-1);
        noiseFreeExcitationNull_same(:,count,:) = noiseFreeExcitationNull_temp;
        noisyExcitationNull_same(:,count,:) = noisyExcitationNull_temp;
    end
end
% Testing Sequence (3 background - 2 test - 3 background)
for f = 1:nFrames
    if f == 4
        [noiseFreeExcitationTest_temp,noisyExcitationTest_temp,~,~,~] = theMosaic.compute(OI_1, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
    elseif f == 5
        [noiseFreeExcitationTest_temp,noisyExcitationTest_temp,~,~,~] = theMosaic.compute(OI_1_2Frames, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
    else
        [noiseFreeExcitationTest_temp,noisyExcitationTest_temp,~,~,~] = theMosaic.compute(OI_2, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
    end
%    [noiseFreeExcitationTest_temp,noisyExcitationTest_temp,~,~,~] = theMosaic.compute(OI_1, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
    for i = 1:nTimebin
        count = i + nTimebin*(f-1);
        noiseFreeExcitationTest_same(:,count,:) = noiseFreeExcitationTest_temp;
        noisyExcitationTest_same(:,count,:) = noisyExcitationTest_temp;
    end
end

%% Photocurrent
noiseFreePhotocurrTest_same = computePhotocurrent(noiseFreeExcitationTest_same, timeAxis8frames, 'none');
noiseFreePhotocurrNull_same = computePhotocurrent(noiseFreeExcitationNull_same, timeAxis8frames, 'none');

noisyPhotocurrTest_same = computePhotocurrent(noisyExcitationTest_same, timeAxis8frames, 'random');
noisyPhotocurrNull_same = computePhotocurrent(noisyExcitationNull_same, timeAxis8frames, 'random');

%% Visualization
visualizeConeResponse(noiseFreeExcitationTest_same, noisyExcitationTest_same, timeAxis8frames, 'excitation');
visualizeAllResponse(noiseFreeExcitationTest_same, noisyExcitationTest_same,noiseFreePhotocurrTest_same, noisyPhotocurrTest_same, timeAxis8frames, 1713, 'Test Stimulus');

visualizeAllResponse(noiseFreeExcitationNull_same, noisyExcitationNull_same,noiseFreePhotocurrNull_same, noisyPhotocurrNull_same, timeAxis8frames, 1713, 'Null Stimulus');

%% SVM
% Simulate a 2AFC task 
taskIntervals = 2;
[classificationMatrix, classLabels] = generateSetUpForClassifier(...
    noisyExcitationTest_same, noisyExcitationNull_same, taskIntervals, 'true');

% Find principal components of the responses
[pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

% Project the responses onto the space formed by the first 4 PC vectors
pcComponentsNumForClassification = 2;
classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);

% Visualize the classification matrix and its projection to the PC space
visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
% Train a binary SVM classifier and visualize the support vectors in 2
% dimensions
svm = fitcsvm(classificationMatrixProjection,classLabels);

% Visualize the data along with the hyperplane computed by the SVM 
visualizeSVMmodel(svm, classificationMatrixProjection, classLabels);
% Measure performance of the SVM classifier using a 10-fold crossvalidation approach
% Perform a 10-fold cross-validation on the trained SVM model
kFold = 10;
CVSVM = crossval(svm,'KFold',kFold);

% Compute classification loss for the in-sample responses using a model 
% trained on out-of-sample responses
fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
% Average percent correct across all folds 
percentCorrect = mean(fractionCorrect)*100


%% Supporting Function
function [meanNoiseFree, meanNoisy, targetCone] = visualizeConeResponse(noiseFreeExcitationResponse, noisyExcitationResponse, timeAxis,responseType)
    figure()
    [~,idx] = max(noiseFreeExcitationResponse(:));
    [~,~,targetConeID] = ind2sub(size(noiseFreeExcitationResponse), idx);
    % Plot the time series response for individual instances
    plot(timeAxis, squeeze(noisyExcitationResponse(:,:,targetConeID)), 'k.');
    hold on;
    % Plot the time series response for the mean of the individual instances
    plot(timeAxis, squeeze(mean(noisyExcitationResponse(:,:,targetConeID),1)), 'g-', 'LineWidth', 2.0);
    % Plot the noise-free time series response in red
    plot(timeAxis, squeeze(mean(noiseFreeExcitationResponse(:,:,targetConeID),1)), 'r', 'LineWidth', 1.5);
    xlabel('time (seconds)');
    if strcmp(responseType,'excitation')
        ylabel('excitations per integration time (R*/cone/tau)');
    elseif strcmp(responseType,'photocurrent')
        ylabel('Photocurrent (pAmps)');
    end
    set(gca, 'FontSize', 16);
    meanNoiseFree = mean(noiseFreeExcitationResponse(:,:,targetConeID));
    meanNoisy = mean(noisyExcitationResponse(:,:,targetConeID));
    targetCone = targetConeID
end

function [meanNoiseFree, meanNoisy] = visualizeSingleConeResponse(noiseFreeExcitationResponse, noisyExcitationResponse, timeAxis, targetConeID, responseType)
    figure()
    % Plot the time series response for individual instances
    plot(timeAxis, squeeze(noisyExcitationResponse(:,:,targetConeID)), '.', 'color', [.5 .5 .5]);
    hold on;
    % Plot the time series response for the mean of the individual instances
    plot(timeAxis, squeeze(mean(noisyExcitationResponse(:,:,targetConeID),1)), 'g-', 'LineWidth', 3.5);
    % Plot the noise-free time series response in red
    plot(timeAxis, squeeze(mean(noiseFreeExcitationResponse(:,:,targetConeID),1)), 'r', 'LineWidth', 2.5);
    xlabel('time (seconds)');
    if strcmp(responseType,'excitation')
        ylabel('excitations per integration time (R*/cone/tau)');
        ylim([40 180])
    elseif strcmp(responseType,'photocurrent')
        ylabel('Photocurrent (pAmps)');
    end
    set(gca, 'FontSize', 16);
    meanNoiseFree = mean(noiseFreeExcitationResponse(:,:,targetConeID));
    meanNoisy = mean(noisyExcitationResponse(:,:,targetConeID));
end

function [meanNoiseFree, meanNoisy] = visualizeAllResponse(noiseFreeExcitationResponse, noisyExcitationResponse, noiseFreePhotocurrResponse, noisyPhotocurrResponse, timeAxis, targetConeID, title_str)
    figure()
    plot(timeAxis, squeeze(noisyExcitationResponse(:,:,targetConeID)), 'c-', 'DisplayName','');
    hold on;
    plot(timeAxis, squeeze(noisyPhotocurrResponse(:,:,targetConeID)), 'm-', 'DisplayName','');
    plot(timeAxis, squeeze(mean(noisyExcitationResponse(:,:,targetConeID),1)), 'b-', 'LineWidth', 2.5, 'DisplayName','Excitation (R*/cone/tau)');
    plot(timeAxis, squeeze(mean(noisyPhotocurrResponse(:,:,targetConeID),1)), 'r-', 'LineWidth', 2.5, 'DisplayName','Photocurrent (pAmps)');
    % plot(timeAxis, squeeze(mean(noiseFreeExcitationResponse(:,:,targetConeID),1)), 'r', 'LineWidth', 1.5);
    xlabel('time (seconds)');
    ylabel('Cone Response');
    ylim([-50 300])
    set(gca, 'FontSize', 16);
    title(title_str)
    meanNoiseFree = mean(noiseFreeExcitationResponse(:,:,targetConeID));
    meanNoisy = mean(noisyExcitationResponse(:,:,targetConeID));
end