%% Generate Tumbling E Object on a Scene
% Two problems unsolved: 
%   1. The function that shifts the tumbling E to specific bar-width
%   2. The function that calcuates the contrast of the E with the background
%       - The background RGB color doesn't seem correct

% Changed from original Tumbling E Scene Engine:
% Use a display Object instead of SPD/ambient .mat file

presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 3);

customSceneParams = struct(...
        'letterRotationDegs', 0, ...            % Letter rotation (0,90,180,270 only)
        'letterHeightPixels', 20, ...           % Letter height in pixels - must be 20
        'letterWidthPixels', 18, ...            % Letter width in pixels - must be 18
        'yPixelsNumMargin', 100, ...             % Y-margin
        'xPixelsNumMargin', 100, ...             % X-margin
        'upSampleFactor', uint8(4), ...         % Upsampling for better centering 
         'chromaSpecification', struct(...
            'type', 'RGBsettings', ...
            'backgroundRGB', [0.4 0.4 0.4], ...   
            'foregroundRGB',  [0.2 0.2 0.2]), ...
        'visualizeScene', false, ...                % Whether to visualize the generated scene
        'displayOBJ', presentationDisplay ...
       );
    %    'spdDataFile', 'BVAMS_White_Guns_At_Max.mat', ...  % Name of file containing the display SPDs
    %    'ambientSPDDataFile', 'BVAMS_White_Background.mat', ... % Name of the file containing the ambient SPD
    %    'plotDisplayCharacteristics', false ...     % Whether to visualize the display characteristics

    % 'chromaSpecification', struct(...
    %               'type', 'chromaLumaLMScontrasts', ...
    %               'backgroundChromaLuma', [0.31 0.32 40], ...
    %               'foregroundLMSConeContrasts', [-0.5 -0.5 0.0]), ...

% text scene can be rotated at 0,90,180,270 degs
deg = 0;
EsceneEngine = createTumblingEsceneEngine(deg,'customSceneParams', customSceneParams);
sceneParamsE = EsceneEngine.sceneComputeFunction();
% Generate scenes with size of 0.5 deg
sizeDegs = 0.5;
EsceneSequence = EsceneEngine.compute(sizeDegs);

dataStruct = sceTumblingEscene(EsceneEngine, sizeDegs, sceneParamsE);
temporalSupport = dataStruct.temporalSupport;

% Get first frame of the scene sequences
Escene = EsceneSequence{1};
% Change mean luminance level
meanLuminanceCdPerM2 = 40;
Escene = sceneAdjustLuminance(Escene, meanLuminanceCdPerM2);

visualizeScene(Escene, ...
            'spatialSupportInDegs', true, ...
            'crossHairsAtOrigin', true, ...
            'displayRadianceMaps', false, ...
            'avoidAutomaticRGBscaling', true, ...
            'noTitle', false);

%% Create Cone Mosaic and initialize params for OI
% 3 frames total = 66 ms (22ms each frame)
integrationTime = 22;
theMosaic = cMosaic('sizeDegs', [1, 1] * sizeDegs, ...
            'eccentricityDegs', [0 0], ...
            'integrationTime', integrationTime / 1000); % [1, 1] * sizeDegs
 
customNeuralResponseParams = struct(...
    'opticsParams',  struct(...
        'type', 'wvf human', ...
        'PolansSubject', 10, ...
        'pupilDiameterMM', 3.0 ...
    ), ...
    'coneMosaicParams', theMosaic); % diaeter must be 1-4 mm
        
%% photopigment excitation
neuralComputeFunction_excitation = @nrePhotopigmentExcitationsCmosaicNoFEM;
theNeuralEngine = neuralResponseEngine(neuralComputeFunction_excitation, customNeuralResponseParams);
instancesNum = 10;
noiseFlags = {'random', 'none'};

[excitationResponses, excitationResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
        EsceneSequence, ...
        temporalSupport, ...
        instancesNum, ...
        'noiseFlags', noiseFlags);
noisyExcitationResponses_0 = excitationResponses('random');
noiseFreeExcitationResponses_0 = excitationResponses('none');

% % Comparison task: 90 deg
% deg = 90;
% EsceneEngine = createTumblingEsceneEngine(deg,'customSceneParams', customSceneParams);
% sceneParamsE = EsceneEngine.sceneComputeFunction();
% EsceneSequence = EsceneEngine.compute(sizeDegs);
% 
% [excitationResponses, excitationResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
%         EsceneSequence, ...
%         temporalSupport, ...
%         instancesNum, ...
%         'noiseFlags', noiseFlags);
% noisyExcitationResponses_90 = excitationResponses('random');
% noiseFreeExcitationResponses_90 = excitationResponses('none');
% 
% % Comparison task: 180 deg
% deg = 180;
% EsceneEngine = createTumblingEsceneEngine(deg,'customSceneParams', customSceneParams);
% sceneParamsE = EsceneEngine.sceneComputeFunction();
% EsceneSequence = EsceneEngine.compute(sizeDegs);
% 
% [excitationResponses, excitationResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
%         EsceneSequence, ...
%         temporalSupport, ...
%         instancesNum, ...
%         'noiseFlags', noiseFlags);
% noisyExcitationResponses_180 = excitationResponses('random');
% noiseFreeExcitationResponses_180 = excitationResponses('none');
% 
% % Comparison task: 270 deg
% deg = 270;
% EsceneEngine = createTumblingEsceneEngine(deg,'customSceneParams', customSceneParams);
% sceneParamsE = EsceneEngine.sceneComputeFunction();
% EsceneSequence = EsceneEngine.compute(sizeDegs);
% 
% [excitationResponses, excitationResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
%         EsceneSequence, ...
%         temporalSupport, ...
%         instancesNum, ...
%         'noiseFlags', noiseFlags);
% noisyExcitationResponses_270 = excitationResponses('random');
% noiseFreeExcitationResponses_270 = excitationResponses('none');


% % merge all frames' responses together
% nTrialsNum = 10;
% nTimebin = 22;
% nFrames = 3;
% ntimeStep = (integrationTime / 1000)/ nTimebin;
% timeAxis3frames = [];
% for t = 0+ntimeStep:ntimeStep:(integrationTime/1000)*nFrames
%     timeAxis3frames(end+1) = t;
% end
% 
% noiseFreeExcitation = zeros(1,nTimebin*nFrames,3822);
% noisyExcitation = zeros(nTrialsNum,nTimebin*nFrames,3822);
% for f = 1:nFrames
%     [excitationResponses, excitationResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
%         EsceneSequence, ...
%         temporalSupport, ...
%         instancesNum, ...
%         'noiseFlags', noiseFlags);
%     noisyExcitationResponses_temp = excitationResponses('random');
%     noiseFreeExcitationResponses_temp = excitationResponses('none');
%     for i = 1:nTimebin
%         count = i + nTimebin*(f-1);
%         noiseFreeExcitation(:,count,:) = noiseFreeExcitationResponses_temp;
%         noisyExcitation(:,count,:) = noisyExcitationResponses_temp;
%     end
% end
%% Visualization
% renderNeuralResponse(theMosaic, noisyExcitationResponses_0);

% renderNeuralResponseSequence(1, noisyExcitation(:,1:22,:), timeAxis3frames(:,1:22), 'noisy');
% renderNeuralResponseSequence(2, noiseFreeExcitation, timeAxis3frames, 'noise-free');
theOI = oiCompute(Escene, theNeuralEngine.neuralPipeline.optics);

[theConeMosaicModulation, theNoisyConeMosaicModulations] = ...
                    theNeuralEngine.neuralPipeline.coneMosaic.compute(theOI, 'nTrials', 4);

domainVisualizationLimits = 0.3*0.5*[-1 1 -1 1];
domainVisualizationTicks = struct('x', -0.2:0.1:0.2, 'y', -0.2:0.1:0.2);
targetWavelengths = [450 550 650];

visualizeRetinalImagesAtWavelengths(theOI, ...
             targetWavelengths, ...
             domainVisualizationLimits, domainVisualizationTicks);

% visualizeRetinalLMSconeImages(theOI, ...
%              domainVisualizationLimits, ...
%              domainVisualizationTicks);

% visualizeConeMosaicActivations(theNeuralEngine.neuralPipeline.coneMosaic, ...
%         theConeMosaicModulation, theNoisyConeMosaicModulations, ...
%          domainVisualizationLimits, ...
%          domainVisualizationTicks)

%% Photocurrent
neuralComputeFunction_photocurrent = @nrePhotocurrentsCmosaicNoFEM;
theNeuralEngine_photocurrent = neuralResponseEngine(neuralComputeFunction_photocurrent, customNeuralResponseParams);
    % It is possible to freeze the noise by specifying a seed for the
    % random number generator through the 'rngSeed' key/value pair.  The
    % compute function should restore the rng to its current state if this
    % is passed, but not otherwise.
[photocurrentResponses, photocurrentResponseTemporalSupportSeconds] = theNeuralEngine_photocurrent.compute(...
        EsceneSequence, ...
        temporalSupport, ...
        instancesNum, ...
        'noiseFlags', noiseFlags); % 'rngSeed', [] ...
noisyPhotocurrentResponse = photocurrentResponses('random');
noiseFreePhotocurrentResponse = photocurrentResponses('none');

% avgFree_excitation = mean(noiseFreeExcitationResponses,"all");
% avgFree_photocurrent = mean(noiseFreePhotocurrentResponse,"all");
% 
% avg_excitation = mean(noisyExcitationResponses,"all");
% avg_photocurrent = mean(noisyPhotocurrentResponse,"all");

%% Test Threshold Engine

% User-supplied computeFunction for the @responseClassifierEngine
classifierComputeFunction = @rcePoissonNWay_OneStimPerTrial;
% classifierComputeFunction = @rcePcaSVMTAFC;

% User-supplied struct with params appropriate for the @responseClassifierEngine computeFunction
customClassifierParams = struct(...
    'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
    'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear 
    'kernelFunction', 'linear', ...     % linear
    'classifierType', 'svm' ...         % binary SVM classifier
    );

% Instantiate a responseClassifierEngine with the above classifierComputeFunctionHandle and custom classifierParams
theClassifierEngine = responseClassifierEngineNWay(classifierComputeFunction, customClassifierParams);

% Train the binary classifier on the above NULL/TEST response set
theResponses = {noisyExcitationResponses_0, noisyExcitationResponses_90,noisyExcitationResponses_180, noisyExcitationResponses_270};

trainingData = theClassifierEngine.compute('train',...
    theResponses);

% Test predictions of the trained classifier on out of sample
% response instances.  This parameter tells us how many out of sample 
% stimulus instances to use.
outOfSampleInstancesNum = 100;

% Repeat out-of-sample predictions a total of N times
    % Compute new respose instances to the NULL and TEST stimuli
%     outOfSampleNullStimResponses = theNeuralEngine.compute(...
%         noisyExcitationResponses_0, ...
%         temporalSupport, ...
%         outOfSampleInstancesNum, ...
%         'noiseFlags', {'random'});
% 
%     outOfSampleTestStimResponses = theNeuralEngine.compute(...
%         noisyExcitationResponses_90, ...
%         temporalSupport, ...
%         outOfSampleInstancesNum, ...
%         'noiseFlags', {'random'});
orientations = [0 90 180 270];
test_data = noisyExcitationResponses_0(1:4,:,:); % ; noisyExcitationResponses_90(1,:,:); noisyExcitationResponses_180(1,:,:); noisyExcitationResponses_270(1,:,:)
% Run the classifier on the new response instances
predictedData = theClassifierEngine.compute('predict',...
    test_data, ...
    orientations);

disp(predictedData.pCorrect);
    % Visualize the classifier performance on the in-sample (training) and 
    % the out of sample responses. Notice that the second component
    % captures almost no variance in the test data. Thats because the principal
    % components are computed on the training data. In that data set the second
    % component is capturing the noise. In the second data set, the
    % second component captures almost zero variance because the noise component
    % is different.
%     debugClassifier = true;
%     if (debugClassifier)
%          plotClassifierResults(trainingData, predictedData);
%          drawnow;
%     end

%% Supporting Functions
function plotClassifierResults(trainingData, predictedData)
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 950 500], 'Color', [1 1 1]);
    minFeature = min([ min(trainingData.features(:)) min(predictedData.features(:)) ]);
    maxFeature = max([ max(trainingData.features(:)) max(predictedData.features(:)) ]);
    
    % The training data
    ax = subplot(1,2,1);
    hold(ax, 'on');
    renderDecisionBoundary(ax,trainingData.decisionBoundary, true); 
    renderFeatures(ax, trainingData.features, trainingData.nominalClassLabels);
    xlabel(ax,'PCA #1 score');
    ylabel(ax,'PCA #2 score');
    title(ax,sprintf('In-sample percent correct: %2.3f', trainingData.pCorrect));
    axis(ax,'square');
    set(ax, 'XLim', [minFeature maxFeature], 'YLim', [minFeature maxFeature], 'FontSize', 12);
    colormap(ax,brewermap(1024, 'RdYlGn'));


    % The predicted data
    ax = subplot(1,2,2);
    hold(ax, 'on');
    renderDecisionBoundary(ax,trainingData.decisionBoundary, false); 
    renderFeatures(ax, predictedData.features, predictedData.nominalClassLabels);
    xlabel(ax,'PCA #1 score');
    ylabel(ax,'PCA #2 score');
    title(ax,sprintf('Out-of-sample percent correct: %2.3f', predictedData.pCorrect));
    axis(ax,'square');
    set(ax, 'XLim', [minFeature maxFeature], 'YLim', [minFeature maxFeature],  'FontSize', 12);
    colormap(ax,brewermap(1024, 'RdYlGn'));
    
end

function renderFeatures(ax, features, nominalLabels)

    idx = find(nominalLabels == 0);
    scatter(ax,features(idx,1), features(idx,2), 64, ...
        'MarkerFaceColor', [0.8 0.8 0.8], ...
        'MarkerEdgeColor', [0.2 0.2 0.2]); 
    hold('on')
    idx = find(nominalLabels == 1);
    scatter(ax,features(idx,1), features(idx,2), 64, ...
        'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', [0.9 0.9 0.9]);

    legend({'decision boundary', 'nominal class 0', 'nomimal class 1'})
end


function renderDecisionBoundary(ax, decisionBoundary, depictStrength)
    if (~isempty(decisionBoundary))
        N = length(decisionBoundary.x);
        if (depictStrength)
        % Decision boundary as a density plot
            imagesc(ax,decisionBoundary.x,decisionBoundary.y,reshape(decisionBoundary.z,[N N]));
        end
        % The decision boundary as a line
        [C,h] = contour(ax,decisionBoundary.x,decisionBoundary.y,reshape(decisionBoundary.z,[N N]), [0 0]);
        h.LineColor = [0 0 0];
        h.LineWidth = 2.0;
    end
end

function renderNeuralResponse(theCMosaic, theResponses)
    hFig = figure(); clf;
    meanResponse = squeeze(mean(theResponses,1));
    ax = subplot(1,2,1);
    theCMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
        'activation', squeeze(theResponses(1,:,:)), 'visualizedConeAperture', 'geometricArea');
    title('1st response instance');
    
    ax = subplot(1,2,2);
    theCMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
        'activation', meanResponse, 'visualizedConeAperture', 'geometricArea');
    title('mean response (%s)');
    colormap(gray);
end

function renderNeuralResponseSequence(figNo, theResponseSequence, theResponseTemporalSupportSeconds, titleLabel)
    figure(figNo); clf;
    meanResponse = squeeze(mean(theResponseSequence,1));
%     cellsNum = size(meanResponse ,1);
    cellsNum = size(theResponseSequence ,1);
    subplot(1,2,1);
    imagesc(theResponseTemporalSupportSeconds, 1:cellsNum, squeeze(theResponseSequence(1,:,:)));
    xlabel('time (sec)');
    ylabel('cells');
    title(sprintf('1st response instance (%s)', titleLabel));
    
    subplot(1,2,2);
    imagesc(theResponseTemporalSupportSeconds, 1:cellsNum, meanResponse);
    xlabel('time (sec)');
    ylabel('cells');
    title(sprintf('mean response (%s)', titleLabel));
    colormap(gray);
end

function visualizeRetinalLMSconeImages(opticalImage, ...
             domainVisualizationLimits, ...
             domainVisualizationTicks)

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 4, ...
       'heightMargin',  0.0, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.00);

    hFig = figure(3);
    clf;
    set(hFig, 'Position', [10 10 1650 500], 'Color', [1 1 1]);

    retinalLMSconeImages = oiGet(opticalImage, 'lms');
    retinalLMSconeImages = retinalLMSconeImages/max(retinalLMSconeImages(:));

    % retrieve the spatial support of the scene(in millimeters)
    spatialSupportMM = oiGet(opticalImage, 'spatial support', 'mm');
    
    % Convert spatial support in degrees
    optics = oiGet(opticalImage, 'optics');
    focalLength = opticsGet(optics, 'focal length');
    mmPerDegree = focalLength*tand(1)*1e3;
    spatialSupportDegs = spatialSupportMM/mmPerDegree;
    spatialSupportX = spatialSupportDegs(1,:,1)*60;
    spatialSupportY = spatialSupportDegs(:,1,2)*60;

    for k = size(retinalLMSconeImages,3):-1:1

        ax = subplot('Position', subplotPosVectors(1,4-k).v);

        imagesc(ax, spatialSupportX, spatialSupportY, squeeze(retinalLMSconeImages(:,:,k)));
        set(ax, 'FontSize', 16);
        axis(ax, 'image'); axis 'xy';
        set(ax, 'XTick', domainVisualizationTicks.x*60, 'YTick', domainVisualizationTicks.y*60, ...
                'XTickLabel', sprintf('%2.2f\n', domainVisualizationTicks.x), ...
                'YTickLabel', sprintf('%2.2f\n', domainVisualizationTicks.y));
        set(ax, 'XLim', domainVisualizationLimits(1:2)*60, 'YLim', domainVisualizationLimits(3:4)*60);
        xlabel('space (deg)');
        
        if (k == 1)
            ylabel('space (deg)');
        end
        colormap(ax, gray(1024));
        hC = colorbar(ax, 'south');
        switch k
            case 1
                hC.Label.String = 'retinal L-cone activation';
            case 2
                hC.Label.String = 'retinal M-cone activation';
            case 3
                hC.Label.String = 'retinal S-cone activation';
        end

        drawnow
    end

%     projectBaseDir = strrep(ISETbioJandJRootPath(), 'toolbox', '');
%     pdfFile = [fullfile(projectBaseDir, 'figures') filesep 'RetinalLMSconeImages.pdf'];
%     NicePlot.exportFigToPDF(pdfFile,hFig, 300);

end

function visualizeRetinalImagesAtWavelengths(opticalImage, ...
             targetWavelengths, ...
             domainVisualizationLimits, ...
             domainVisualizationTicks)

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 4, ...
       'heightMargin',  0.0, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.00);

    wavelengthSupport = oiGet(opticalImage, 'wave');
    hFig = figure(3);
    clf;
    set(hFig, 'Position', [10 10 1650 500], 'Color', [1 1 1]);

    retinalImagePhotonRate = oiGet(opticalImage, 'photons');

    % retrieve the spatial support of the scene(in millimeters)
    spatialSupportMM = oiGet(opticalImage, 'spatial support', 'mm');
    
    % Convert spatial support in degrees
    optics = oiGet(opticalImage, 'optics');
    focalLength = opticsGet(optics, 'focal length');
    mmPerDegree = focalLength*tand(1)*1e3;
    spatialSupportDegs = spatialSupportMM/mmPerDegree;
    spatialSupportX = spatialSupportDegs(1,:,1)*60;
    spatialSupportY = spatialSupportDegs(:,1,2)*60;

    for k = 1:numel(targetWavelengths)
        [~,wIndex] = min(abs(wavelengthSupport-targetWavelengths(k)));
        targetWavelength = wavelengthSupport(wIndex);
    
        ax = subplot('Position', subplotPosVectors(1,k).v);
        imagesc(ax, spatialSupportX, spatialSupportY, squeeze(retinalImagePhotonRate(:,:,wIndex)));
        set(ax, 'FontSize', 16);
        axis(ax, 'image'); axis 'xy';
        set(ax, 'XTick', domainVisualizationTicks.x*60, 'YTick', domainVisualizationTicks.y*60, ...
                'XTickLabel', sprintf('%2.2f\n', domainVisualizationTicks.x), ...
                'YTickLabel', sprintf('%2.2f\n', domainVisualizationTicks.y));
        set(ax, 'XLim', domainVisualizationLimits(1:2)*60, 'YLim', domainVisualizationLimits(3:4)*60);
        xlabel('space (deg)');
        
        if (k == 1)
            ylabel('space (deg)');
        end
        colormap(ax, gray(1024));
        hC = colorbar(ax, 'south');
        hC.Label.String = sprintf('retinal irradiance (photons/m^2/sec)');

        title(ax, sprintf('%d nm', targetWavelength));
        drawnow
    end

%     projectBaseDir = strrep(ISETbioJandJRootPath(), 'toolbox', '');
%     pdfFile = [fullfile(projectBaseDir, 'figures') filesep 'RetinalImages.pdf'];
%     NicePlot.exportFigToPDF(pdfFile,hFig, 300);

end

function visualizeConeMosaicActivations(theConeMosaic, ...
        theConeMosaicModulation, theNoisyConeMosaicModulations, ...
         domainVisualizationLimits, ...
         domainVisualizationTicks)

     domainVisualizationTicks.y = [-0.1 0 0.1];
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 4, ...
       'heightMargin',  0.0, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.00);

    hFig = figure(4);
    clf;
    set(hFig, 'Position', [10 10 1650 500], 'Color', [1 1 1]);


    for iOri = 1:4
        ax = subplot('Position', subplotPosVectors(1,iOri).v);
        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'activation', theConeMosaicModulation{iOri}, ...
            'activationRange', 20*[-1 1], ...
            'verticalActivationColorBarInside', true, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'noYLabel', (iOri>1), ...
            'plotTitle', sprintf('cone modulation (max: %2.1f%%)',max(abs(theConeMosaicModulation{iOri}(:)))), ...
            'fontSize', 16, ...
            'backgroundColor', [0 0 0]);
        xtickangle(ax, 0);
        ytickangle(ax, 0);
    end

    projectBaseDir = strrep(ISETbioJandJRootPath(), 'toolbox', '');
    pdfFile = [fullfile(projectBaseDir, 'figures') filesep 'MeanConeMosaicActivations.pdf'];
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);


    hFig = figure(5);
    clf;
    set(hFig, 'Position', [10 10 1650 500], 'Color', [1 1 1]);

    for iInstance = 1:4
        ax = subplot('Position', subplotPosVectors(1,iInstance).v);
        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'activation', squeeze(theNoisyConeMosaicModulations{1}(iInstance,:,:)), ...
            'activationRange', 20*[-1 1], ...
            'verticalActivationColorBarInside', true, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'noYLabel', (iOri>1), ...
            'plotTitle', sprintf('noisy cone mosaic activation\n(instance %d)', iInstance),...
            'fontSize', 16, ...
            'backgroundColor', [0 0 0]);
        xtickangle(ax, 0);
        ytickangle(ax, 0);
    end

    projectBaseDir = strrep(ISETbioJandJRootPath(), 'toolbox', '');
    pdfFile = [fullfile(projectBaseDir, 'figures') filesep 'NoisyConeMosaicActivations.pdf'];
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);


    hFig = figure(6);
    clf;
    set(hFig, 'Position', [10 10 1650 500], 'Color', [1 1 1]);

    ax = subplot('Position', subplotPosVectors(1,4).v);
    theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'plotTitle', 'cone mosaic',...
            'fontSize', 16, ...
            'backgroundColor', [0 0 0]);

    wave = theConeMosaic.wave;
    photopigment = theConeMosaic.pigment;

    ax = subplot('Position', subplotPosVectors(1,1).v);
    plot(ax, wave, photopigment.quantalEfficiency(:,3), 'bo-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
    xlabel(ax, 'wavelength (nm)');
    ylabel(ax, 'quantal efficiency');
    axis(ax, 'square'); grid(ax, 'on');
    set(ax, 'YLim', [0 0.5], 'YTick', 0:0.1:0.5, 'FontSize', 16)
    title(ax, 'S-cone');

    ax = subplot('Position', subplotPosVectors(1,2).v);
    plot(ax, wave, photopigment.quantalEfficiency(:,2), 'go-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [0.5 1 0.5], 'MarkerEdgeColor', [0.5 0 0], 'LineWidth', 1.5);
    xlabel(ax, 'wavelength (nm)');
    axis(ax, 'square'); grid(ax, 'on');
    set(ax, 'YLim', [0 0.5], 'YTick', 0:0.1:0.5, 'FontSize', 16)
    title(ax, 'M-cone');

    ax = subplot('Position', subplotPosVectors(1,3).v);
    plot(ax, wave, photopigment.quantalEfficiency(:,1), 'ro-', ...
        'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5],  'LineWidth', 1.5);
    xlabel(ax, 'wavelength (nm)');
    axis(ax, 'square'); grid(ax, 'on');
    set(ax, 'YLim', [0 0.5], 'YTick', 0:0.1:0.5, 'FontSize', 16)
    title(ax, 'L-cone');

    projectBaseDir = strrep(ISETbioJandJRootPath(), 'toolbox', '');
    pdfFile = [fullfile(projectBaseDir, 'figures') filesep 'ConeMosaic.pdf'];
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);

end
   

    
  

        
        
        
        
        
        
        