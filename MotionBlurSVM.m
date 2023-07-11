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
    'meanLuminanceCdPerM2', 20, ...
    'center', [0 0], ...
    'pixelsAlongWidthDim', [], ...          
    'pixelsAlongHeightDim', [] ...          
    );
scene1 = generateGaborSceneCustomize(...
    'stimParams', stimParams,...
    'presentationDisplay', presentationDisplay, ...
    'minimumPixelsPerHalfPeriod', 10);

stimParams.center = [0 20];
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
integrationTime = 50;
theMosaic = cMosaic('sizeDegs', [1, 1] * stimParams.sizeDegs, ...
            'integrationTime', integrationTime / 1000);

%% 8 Frames Experiment: 3 - 2 - 3
nTrialsNum = 10;
nTimebin = 50;
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
% visualizeConeResponse(noiseFreeExcitationNull_same, noisyExcitationNull_same, timeAxis8frames, 'excitation');
% visualizeSingleConeResponse(noiseFreeExcitationTest_same, noisyExcitationTest_same, timeAxis8frames, 1713, 'excitation');
% visualizeSingleConeResponse(noiseFreeExcitationNull_same, noisyExcitationNull_same, timeAxis8frames, 1713, 'excitation');
% Visualize Photocurrent
% visualizeSingleConeResponse(noiseFreePhotocurrTest_same, noisyPhotocurrTest_same, timeAxis8frames, 1713, 'photocurrent');
% visualizeSingleConeResponse(noiseFreePhotocurrNull_same, noisyPhotocurrNull_same, timeAxis8frames, 1713, 'photocurrent');
% Try All
visualizeAllResponse(noiseFreeExcitationTest_same, noisyExcitationTest_same,noiseFreePhotocurrTest_same, noisyPhotocurrTest_same, timeAxis8frames, 1713, 'Test Stimulus');

visualizeAllResponse(noiseFreeExcitationNull_same, noisyExcitationNull_same,noiseFreePhotocurrNull_same, noisyPhotocurrNull_same, timeAxis8frames, 1713, 'Null Stimulus');

%% Compute Cone Excitation
nTrialsNum = 20;
% duration = 100/1000;
% without eye movement, 2 frames
% theMosaic.emGenSequence(duration, 'nTrials', 1, 'centerPaths', true,'microsaccadeType', 'none');
nTimebin = 5;
[noiseFreeExcitationStimulus_c1_f1,noiseFreeExcitationStimulus_c1_f2,...
    noiseFreeExcitationNull_c1_f1,noiseFreeExcitationNull_c1_f2] = deal(zeros(1,nTimebin,3822));
[noisyExcitationStimulus_c1_f1, noisyExcitationStimulus_c1_f2, ...
    noisyExcitationNull_c1_f1, noisyExcitationNull_c1_f2] = deal(zeros(nTrialsNum,nTimebin,3822));

[noiseFreeExcitationInstanceStimulus_c1_f1,noisyExcitationInstanceStimulus_c1_f1,~,~,~] = theMosaic.compute(OI_1, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
[noiseFreeExcitationInstanceNull_c1_f1,noisyExcitationInstanceNull_c1_f1,~,~,~] = theMosaic.compute(OI_2, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);

[noiseFreeExcitationInstanceStimulus_c1_f2,noisyExcitationInstanceStimulus_c1_f2,~,~,~] = theMosaic.compute(OI_1, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
[noiseFreeExcitationInstanceNull_c1_f2,noisyExcitationInstanceNull_c1_f2,~,~,~] = theMosaic.compute(OI_2, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);

for i = 1:nTimebin
    noiseFreeExcitationStimulus_c1_f1(:,i,:) = noiseFreeExcitationInstanceStimulus_c1_f1;
    noiseFreeExcitationNull_c1_f1(:,i,:) = noiseFreeExcitationInstanceNull_c1_f1;
    noiseFreeExcitationStimulus_c1_f2(:,i,:) = noiseFreeExcitationInstanceStimulus_c1_f2;
    noiseFreeExcitationNull_c1_f2(:,i,:) = noiseFreeExcitationInstanceNull_c1_f2;
    noisyExcitationStimulus_c1_f1(:,i,:) = noisyExcitationInstanceStimulus_c1_f1;
    noisyExcitationNull_c1_f1(:,i,:) = noisyExcitationInstanceNull_c1_f1;
    noisyExcitationStimulus_c1_f2(:,i,:) = noisyExcitationInstanceStimulus_c1_f2;
    noisyExcitationNull_c1_f2(:,i,:) = noisyExcitationInstanceNull_c1_f2;
end
noiseFreeExcitationStimulus_c1 = [noiseFreeExcitationStimulus_c1_f1 noiseFreeExcitationStimulus_c1_f2];
noisyExcitationStimulus_c1 = [noisyExcitationStimulus_c1_f1 noisyExcitationStimulus_c1_f2];
noiseFreeExcitationNull_c1 = [noiseFreeExcitationNull_c1_f1 noiseFreeExcitationNull_c1_f2];
noisyExcitationNull_c1 = [noisyExcitationNull_c1_f1 noisyExcitationNull_c1_f2];

ntimeStep = (integrationTime / 1000)/ nTimebin;
timeAxis = [];
for t = 0+ntimeStep:ntimeStep:(integrationTime/1000)
    timeAxis(end+1) = t;
end
timeAxis2frames = [];
for t = 0+ntimeStep:ntimeStep:(integrationTime/1000)*2
    timeAxis2frames(end+1) = t;
end

% Average Exitations for comparison
avgFreeStimulus_c1_f1 = mean(noiseFreeExcitationInstanceStimulus_c1_f1);
avgFreeNull_c1_f1 = mean(noiseFreeExcitationInstanceNull_c1_f1);
avgFreeStimulus_c1_f2 = mean(noiseFreeExcitationInstanceStimulus_c1_f2);
avgFreeNull_c1_f2 = mean(noiseFreeExcitationInstanceNull_c1_f2);

avgStimulus_c1_f1 = mean(noisyExcitationInstanceStimulus_c1_f1,"all");
avgNull_c1_f1 = mean(noisyExcitationInstanceNull_c1_f1,"all");
avgStimulus_c1_f2 = mean(noisyExcitationInstanceStimulus_c1_f2,"all");
avgNull_c1_f2 = mean(noisyExcitationInstanceNull_c1_f2,"all");

figure
x = [1 2];
y1 = [avgFreeStimulus_c1_f1 avgFreeStimulus_c1_f2];
y2 = [avgFreeNull_c1_f1 avgFreeNull_c1_f2];
y3 = [avgStimulus_c1_f1 avgStimulus_c1_f2];
y4 = [avgNull_c1_f1 avgNull_c1_f2];
box_stimulus = [avgFreeStimulus_c1_f1 avgFreeStimulus_c1_f2 avgStimulus_c1_f1 avgStimulus_c1_f2];
box_null = [avgFreeNull_c1_f1 avgFreeNull_c1_f2 avgNull_c1_f1 avgNull_c1_f2];
box_free = [avgFreeStimulus_c1_f1 avgFreeNull_c1_f1];

group = [ones(size(box_stimulus)) 2 * ones(size(box_null))];
h = boxplot([box_stimulus box_null], group, 'Color','b');
set(gca,'XTickLabel',{'Test','Null'},'linew',1)
set(h,'LineWidth',2);
% boxchart([1 2],box_data)
hold on;
plot(x, box_free, '.r', 'markersize', 20);
% % boxchart([2],box_null)
% % plot(x, y1, 'or');
% % plot(x, y2, 'ob');
% % plot(x, y3, 'or');
% % plot(x, y4, 'ob');
% % xlim([0.5 2.5])
% % xticks([1 2])
% title('Cone Excitations (no FEM)','fontsize', 15)
% xlabel('Frames')
ylabel('Cone Excitation', 'fontsize', 12)
% legend('stimulus','null')
% % legend('noise free stimulus','noise free null','noisy stimulus','noisy null')
hold off

% visualizeConeExcitation(noiseFreeExcitationStimulus_c1, noisyExcitationStimulus_c1, timeAxis2frames, 'true');
% visualizeConeExcitation(noiseFreeExcitationNull_c1, noisyExcitationNull_c1, timeAxis2frames, 'true');

% theMosaic.visualize('activation',noiseFreeExcitationInstanceStimulus_c1_f1,'plotTitle','Test Stimulus','verticalActivationColorBar', true);
% theMosaic.visualize('activation',noiseFreeExcitationInstanceNull_c1_f1,'plotTitle','Null Stimulus','verticalActivationColorBar', true);
% Compute Photocurrent Try
keptDuration = 50/1000;
noiseFreePhotocurrStimulus_c1_f1 = computePhotocurrent(noiseFreeExcitationStimulus_c1_f1, timeAxis, 'none', keptDuration);
noiseFreePhotocurrStimulus_c1_f2 = computePhotocurrent(noiseFreeExcitationStimulus_c1_f2, timeAxis, 'none', keptDuration);
noiseFreePhotocurrNull_c1_f1 = computePhotocurrent(noiseFreeExcitationNull_c1_f1, timeAxis, 'none', keptDuration);
noiseFreePhotocurrNull_c1_f2 = computePhotocurrent(noiseFreeExcitationNull_c1_f2, timeAxis, 'none', keptDuration);

noisyPhotocurrStimulus_c1_f1 = computePhotocurrent(noisyExcitationStimulus_c1_f1, timeAxis, 'random', keptDuration);
noisyPhotocurrStimulus_c1_f2 = computePhotocurrent(noisyExcitationStimulus_c1_f2, timeAxis, 'random', keptDuration);
noisyPhotocurrNull_c1_f1 = computePhotocurrent(noisyExcitationNull_c1_f1, timeAxis, 'random', keptDuration);
noisyPhotocurrNull_c1_f2 = computePhotocurrent(noisyExcitationNull_c1_f2, timeAxis, 'random', keptDuration);

avgPhotocurrFreeStimulus_c1_f1 = mean(noiseFreePhotocurrStimulus_c1_f1,"all");
avgPhotocurrFreeNull_c1_f1 = mean(noiseFreePhotocurrNull_c1_f1,"all");
avgPhotocurrFreeStimulus_c1_f2 = mean(noiseFreePhotocurrStimulus_c1_f2,"all");
avgPhotocurrFreeNull_c1_f2 = mean(noiseFreePhotocurrNull_c1_f2,"all");

avgPhotocurrStimulus_c1_f1 = mean(noisyPhotocurrStimulus_c1_f1,"all");
avgPhotocurrNull_c1_f1 = mean(noisyPhotocurrNull_c1_f1,"all");
avgPhotocurrStimulus_c1_f2 = mean(noisyPhotocurrStimulus_c1_f2,"all");
avgPhotocurrNull_c1_f2 = mean(noisyPhotocurrNull_c1_f2,"all");

% figure
% x = [1 2];
% y1 = [avgPhotocurrFreeStimulus_c1_f1 avgPhotocurrFreeStimulus_c1_f2];
% y2 = [avgPhotocurrFreeNull_c1_f1 avgPhotocurrFreeNull_c1_f2];
% y3 = [avgPhotocurrStimulus_c1_f1 avgPhotocurrStimulus_c1_f2];
% y4 = [avgPhotocurrNull_c1_f1 avgPhotocurrNull_c1_f2];
% plot(x, y1, 'or');
% hold on;
% plot(x, y2, 'ob');
% plot(x, y3, 'or');
% plot(x, y4, 'ob');
% xlim([0.5 2.5])
% xticks([1 2])
% title('Photocurrent (no FEM)','fontsize', 15)
% xlabel('Frames')
% ylabel('Cone Excitation')
% legend('noise free stimulus','noise free null','noisy stimulus','noisy null')
% hold off

figure
x = [1 2];
box_stimulus = [avgPhotocurrFreeStimulus_c1_f1 avgPhotocurrFreeStimulus_c1_f2 avgPhotocurrStimulus_c1_f1 avgPhotocurrStimulus_c1_f2];
box_null = [avgPhotocurrFreeNull_c1_f1 avgPhotocurrFreeNull_c1_f2 avgPhotocurrNull_c1_f1 avgPhotocurrNull_c1_f2];
box_free = [avgPhotocurrFreeStimulus_c1_f1 avgPhotocurrFreeNull_c1_f1];

group = [ones(size(box_stimulus)) 2 * ones(size(box_null))];
h = boxplot([box_stimulus box_null], group, 'Color','b');
set(gca,'XTickLabel',{'Test','Null'},'linew',1)
set(h,'LineWidth',2);
hold on;
plot(x, box_free, '.r', 'markersize', 20);
ylabel('Photocurrent (pAmps)', 'fontsize', 12)
hold off

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
    ylim([-100 1200])
    set(gca, 'FontSize', 16);
    title(title_str)
    meanNoiseFree = mean(noiseFreeExcitationResponse(:,:,targetConeID));
    meanNoisy = mean(noisyExcitationResponse(:,:,targetConeID));
end