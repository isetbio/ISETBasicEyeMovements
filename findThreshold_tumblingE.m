
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 1.5);
params = struct(...
    'letterSizesNumExamined', 15, ...                           % How many sizes to use for sampling the psychometric curve
    'maxLetterSizeDegs', 0.5, ...                               % The maximum letter size in degrees of visual angle
    'sceneUpSampleFactor', uint8(4), ...                        % Upsample scene, so that the pixel for the smallest scene is < cone aperture
    'mosaicIntegrationTimeSeconds', 22/1000, ...                % Integration time, here 300 msec
    'nTest', 512, ...                                           % Number of trial to use for computing Pcorrect, 512
    'thresholdP', 0.80, ...                                     % Probability correct level for estimating threshold performance
    'chromaSpecification', struct(...
            'type', 'RGBsettings', ...
            'backgroundRGB', [0.4 0.4 0.4], ...   
            'foregroundRGB',  [0.2 0.2 0.2]), ...
    'visualizedPSFwavelengths', [], ... % 380:10:770, ...        % Vector with wavelengths for visualizing the PSF. If set to empty[] there is no visualization.
    'visualizeDisplayCharacteristics', ~true, ...               % Flag, indicating whether to visualize the display characteristics
    'visualizeScene', ~true, ...                                 % Flag, indicating whether to visualize one of the scenes
    'displayOBJ', presentationDisplay ...
);
% 'spdDataFile', 'BVAMS_White_Guns_At_Max.mat', ...           % Datafile containing the display SPDs
% 'psfDataSubDir', 'FullVis_PSFs_10nm_Subject9', ...          % Subdir where the PSF data live
% 'psfDataFile', '',...                                       % Datafile containing the PSF data
    
% examinedPSFDataFiles = {...
%     'Uniform_FullVis_LCA_zero_TCA_zero.mat' ...
%     'Uniform_FullVis_LCA_low_TCA_zero.mat' ...
%     'Uniform_FullVis_LCA_high_TCA_zero.mat' ...
%     'Uniform_FullVis_LCA_zero_TCA_low.mat' ...
%     'Uniform_FullVis_LCA_low_TCA_low.mat' ...
%     'Uniform_FullVis_LCA_high_TCA_low.mat' ...
%     'Uniform_FullVis_LCA_zero_TCA_high.mat' ...
%     'Uniform_FullVis_LCA_low_TCA_high.mat' ...
%     'Uniform_FullVis_LCA_high_TCA_high.mat' ...
%     };


theConeMosaic = [];
for iPSF = 1:1% numel(examinedPSFDataFiles{1})
    tic
%     params.psfDataFile = examinedPSFDataFiles{iPSF};
    theConeMosaic = runSimulation(params, theConeMosaic);
    toc
end

function theConeMosaic = runSimulation(params, theConeMosaic)
    meanLuminance = 5;
    % Unpack simulation params
    letterSizesNumExamined = params.letterSizesNumExamined;
    maxLetterSizeDegs = params.maxLetterSizeDegs;
    mosaicIntegrationTimeSeconds = params.mosaicIntegrationTimeSeconds;
    nTest = params.nTest;
    thresholdP = params.thresholdP;
%     spdDataFile = params.spdDataFile;
%     psfDataFile = fullfile(params.psfDataSubDir, params.psfDataFile);

    % Visualization of the PSF stack
%     if (~isempty(params.visualizedPSFwavelengths))
%         psfRangeArcMin = 10;
%         visualizePSFstack(theCustomPSFOptics, params.visualizedPSFwavelengths, psfRangeArcMin)
%     end
    
    if (isempty(theConeMosaic))
        mosaicSizeDegs = maxLetterSizeDegs*1*[1 1]; % 1.25
        theConeMosaic = cMosaic('sizeDegs', mosaicSizeDegs, ...
                'eccentricityDegs', [0 0], ...
                'integrationTime', mosaicIntegrationTimeSeconds); 
    end 

    %% Create neural response engine
    % TO DO: change it to original Neural Response Engine
    customNeuralResponseParams = struct(...
        'opticsParams',  struct(...
            'type', 'wvf human', ...
            'PolansSubject', 10, ...
            'pupilDiameterMM', 3.0 ...
        ), ...
        'coneMosaicParams', theConeMosaic); % diameter must be 1-4 mm   
    
    neuralComputeFunction_excitation = @nrePhotopigmentExcitationsCmosaicNoFEM; % @nrePhotopigmentExcitationsCmosaicSingleShot
    theNeuralEngine = neuralResponseEngine(neuralComputeFunction_excitation, customNeuralResponseParams);
%     noiseFlags = {'random', 'none'};
%     theOI = oiCompute(Escene, theNeuralEngine.neuralPipeline.optics);
    
%     theNeuralEngine.customNeuralPipeline(struct(...
%               'coneMosaic', theConeMosaic, ...
%               'optics', theCustomPSFOptics)); % TO DO: change coneMosaic & optics (check: coneMosaic/cMosaic)

    % Poisson n-way AFC
    classifierEngine = responseClassifierEngineNWay(@rcePoissonNWay_OneStimPerTrial);
    % Parameters associated with use of the Poisson classifier.
    classifierPara = struct('trainFlag', 'none', ...
                            'testFlag', 'random', ...
                            'nTrain', 1, 'nTest', nTest);

    % Tumbling E setup
    orientations = [0 90 180 270];
    
    %% Parameters for threshold estimation/quest engine
    thresholdParameters = struct(...
        'maxParamValue', maxLetterSizeDegs, ...    % The maximum value of the examined param (letter size in degs)
        'logThreshLimitLow', 2.0, ...              % minimum log10(normalized param value)
        'logThreshLimitHigh', 0.0, ...             % maximum log10(normalized param value)
        'logThreshLimitDelta', 0.01, ...
        'slopeRangeLow', 1/20, ...
        'slopeRangeHigh', 500/20, ...
        'slopeDelta', 2/20, ...
        'thresholdCriterion', thresholdP, ...
        'guessRate', 1/numel(orientations), ...
        'lapseRate', [0 0.02]);
    
    % Parameters for Quest
    questEnginePara = struct( ...
        'qpPF',@qpPFWeibullLog, ...
        'minTrial', nTest*letterSizesNumExamined, ...
        'maxTrial', nTest*letterSizesNumExamined, ...
        'numEstimator', 1, ...
        'stopCriterion', 0.05, ...
        'employMethodOfConstantStimuli', ~true, ...
        'nTest', 64);

    % Generate a customSceneParams struct from the defaultSceneParams
    % so we can set certain scene params of interest
    theSceneEngine = createTumblingEsceneEngine(0);
    customSceneParams = theSceneEngine.sceneComputeFunction();
    customSceneParams.yPixelsNumMargin = 100;
    customSceneParams.xPixelsNumMargin = 100;
    customSceneParams.chromaSpecification = params.chromaSpecification;
    customSceneParams.displayOBJ = params.displayOBJ;
    customSceneParams.upSampleFactor = params.sceneUpSampleFactor;
        
    if (params.visualizeDisplayCharacteristics)
        visualizationSceneParams = customSceneParams;
        visualizationSceneParams.plotDisplayCharacteristics = true;
        theSceneEngine.sceneComputeFunction(theSceneEngine,0.01, visualizationSceneParams);
    end
    
    if (params.visualizeScene)
        visualizationSceneParams = customSceneParams;
        visualizationSceneParams.visualizeScene = true;
        theSceneEngine.sceneComputeFunction(theSceneEngine,0.01, visualizationSceneParams);
    end
    
     clear 'theSceneEngine';
    
    % Generate scene engines for the tumbling E's (4 orientations)
    tumblingEsceneEngines = cell(1,numel(orientations));
    for iOri = 1:numel(orientations)
       tumblingEsceneEngines{iOri} = createTumblingEsceneEngine(...
           orientations(iOri), ...
           'customSceneParams', customSceneParams);
    end

    % Generate background scene engine.params for the background scene
    sceneParams = tumblingEsceneEngines{1}.sceneComputeFunction();
    backgroundSceneParams = sceneParams;
    backgroundSceneParams.chromaSpecification.foregroundRGB = sceneParams.chromaSpecification.backgroundRGB;
    backgroundSceneEngine = createTumblingEsceneEngine(orientations(1), 'customSceneParams', backgroundSceneParams);

    % Compute psychometric function for the 4AFC paradigm with the 4 E scenes
    [threshold, questObj, psychometricFunction, fittedPsychometricParams] = computeParameterThreshold(...
            tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
            classifierPara, thresholdParameters, questEnginePara, meanLuminance, ...
            'visualizeAllComponents', ~true, ...
            'beVerbose', true);

    % Plot the derived psychometric function
%     pdfFileName = sprintf('Performance_%s_Reps_%d.pdf', strrep(params.psfDataFile, '.mat', ''), nTest);
    plotDerivedPsychometricFunction(questObj, threshold, fittedPsychometricParams, ...
        thresholdParameters); % , 'xRange', [0.02 0.2]
    
    disp('curve already displayed')

%     pdfFileName = sprintf('Simulation_%s_Reps_%d.pdf', strrep(params.psfDataFile, '.mat', ''), nTest);
%     visualizeSimulationResults(questObj, threshold, fittedPsychometricParams, ...
%         thresholdParameters, tumblingEsceneEngines, theNeuralEngine);

    % Export the results
%     exportFileName = sprintf('Results_%s_Reps_%d.mat', strrep(params.psfDataFile, '.mat', ''), nTest);
%     exportSimulation(questObj, threshold, fittedPsychometricParams, ...
%         thresholdParameters, classifierPara, questEnginePara, ...
%         tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
%         exportFileName);

end
