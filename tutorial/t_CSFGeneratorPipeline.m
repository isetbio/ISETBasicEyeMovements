%% CSF Generator Try
% Stimulus 1 
sceneComputeFunction = @sceGrating; % @sceUniformFieldTemporalModulation
sceneParams = struct(...
    'viewingDistanceMeters', 0.5, ...              % display: viewing distance
    'bitDepth', 20, ...                             % display: length of LUT
    'gammaTableExponent', 2.0, ...                  % display: shape of LUT, 1 = Linear
    'spectralSupport', 400:10:750, ...              % display: spectral support of the primary SPDs, in nanometers
    'meanLuminanceCdPerM2', 15, ...                 % background: mean luminance, in candellas per meter squared
    'meanChromaticityXY', [0.3 0.32], ...           % background: neutral mean chromaticity
    'coneContrastModulation', [0.09 -0.09 0.0], ... % chromatic direction: LMS cone contrasts
    'warningInsteadOfErrorOnOutOfGamut', false, ... % chromatic: whether to throw an error or a warning if out-of-gamut pixels are generated
    'fovDegs', 0.5, ...                             % spatial: field of view, in degrees
    'minPixelsNumPerCycle', 10, ...                 % spatial: min number of pixels per grating spatial period (to minimize aliasing effects)
    'pixelsNum', 256, ...                           % spatial: desired size of stimulus in pixels - this could be larger depending on the minPixelsNumPerCycle, the SF, and the FOV
    'spatialFrequencyCyclesPerDeg', 10, ...          % spatial: grating spatial freequency, in cycles/deg
    'orientationDegs', 0, ...                       % spatial: grating orientation, in degrees
    'spatialPhaseDegs', 180, ...                     % spatial: grating spatial phase, in degrees
    'spatialPositionDegs', [0.0 0.0], ...           % spatial: center of grating, in degrees
    'spatialModulation', 'harmonic', ...            % spatial: contrast modulation - choose between {'harmonic', 'square'} for sinusoidal/square spatial modulation 
    'spatialModulationDomain', 'cartesian', ...     % spatial: domain of spatial modulation - choose between {'cartesian', 'polar'}
    'spatialEnvelope', 'soft', ...                  % spatial: envelope - choose between {'disk', 'rect', 'soft'}
    'spatialEnvelopeRadiusDegs', 0.2/3, ...          % spatial: radius of the spatial envelope, in degs    
    'temporalModulation', 'drifted', ...            % temporal modulation mode: choose between {'flashed', 'drifted', 'counter phase modulated'}
    'temporalModulationParams', struct(...          % temporal: modulation params struct
        'temporalFrequencyHz', 50, ...            %   params relevant to the temporalModulationMode
        'stimDurationTemporalCycles', 5), ...            %   params relevant to the temporalModulationMode
    'frameDurationSeconds', 20/1000 ...            % temporal: frame duration, in seconds
    );
testContrast = 0.08;
nullContrast = 0.0;
theSceneEngine = sceneEngine(sceneComputeFunction, sceneParams);
[theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
% theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
% theSceneEngine.visualizeStaticFrame(theSceneSequence); % ,'frameToVisualize', 1

% create cone mosaic
integrationTime = 20;
theMosaic = cMosaic('sizeDegs', [1, 1] * sceneParams.fovDegs, ...
            'eccentricityDegs', [0 0], ...
            'integrationTime', integrationTime / 1000);
        
customNeuralResponseParams = struct(...
    'opticsParams',  struct(...
        'type', 'wvf human', ...
        'PolansSubject', 10, ...
        'pupilDiameterMM', 3.0 ...
    ), ...
    'eyeMovementsParams', struct(...
        'driftModel', 'default', ...
        'durationSeconds', 300/1000, ...
        'keptResponsesDurationSeconds', 6000/1000 ...
    ), ...
    'coneMosaicParams', theMosaic);

% cone excitation
% The neural compute function to employ
neuralComputeFunction_excitation = @nrePhotopigmentExcitationsCmosaicEyeMovements;
theNeuralEngine = neuralResponseEngine(neuralComputeFunction_excitation, customNeuralResponseParams);
instancesNum = 1;
noiseFlags = {'random', 'none'};
[excitationResponses, excitationResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
        theSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        instancesNum, ...
        'noiseFlags', noiseFlags);
noisyExcitationResponses = excitationResponses('random');
noiseFreeExcitationResponses = excitationResponses('none');

% Instantiate a neuralResponseEngine with the custom neural response params
neuralComputeFunction_photocurrent = @nrePhotocurrentsCmosaicEyeMovements;
theNeuralEngine_photocurrent = neuralResponseEngine(neuralComputeFunction_photocurrent, customNeuralResponseParams);
% Compute instances of neural responses to the input scene sequence.
% It is possible to freeze the noise by specifying a seed for the
% randome number generator through the 'rngSeed' key/value pair.  The
% compute function should restore the rng to its current state if this
% is passed, but not otherwise.
[photocurrentResponses, photocurrentResponseTemporalSupportSeconds] = theNeuralEngine_photocurrent.compute(...
        theSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        instancesNum, ...
        'noiseFlags', noiseFlags); % 'rngSeed', [] ...
noisyPhotocurrentResponse = photocurrentResponses('random');
noiseFreePhotocurrentResponse = photocurrentResponses('none');

avgFree_excitation = mean(noiseFreeExcitationResponses,"all");
avgFree_photocurrent = mean(noiseFreePhotocurrentResponse,"all");

avg_excitation = mean(noisyExcitationResponses,"all");
avg_photocurrent = mean(noisyPhotocurrentResponse,"all");
