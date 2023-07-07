%% Create a display object
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 0.50);

%% Create a scene object for stimulus

% specify stimulus parameters
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 12, ... % 12 cycles/deg
    'orientationDegs', 90, ...              % 90 degrees
    'phaseDegs', 90, ...                    % spatial phase degrees
    'sizeDegs', 0.5, ...                    % 0.5 x 0.5 size
    'sigmaDegs', 0.2/3, ...                 % sigma of Gaussian envelope, in degrees
    'contrast', 0.8,...                     % 0.8 Michelson contrast
    'meanLuminanceCdPerM2', 40, ...         % 40 cd/m2 mean luminance
    'pixelsAlongWidthDim', [], ...          % pixels- width dimension
    'pixelsAlongHeightDim', [] ...          % pixel- height dimension
    );
% compute the pixels that the stimulus will occupy

% retrieve the display's pixel size 
displayPixelSizeMeters = displayGet(presentationDisplay, 'sample spacing');
% retrieve the display's viewing distance
viewingDistanceMeters = displayGet(presentationDisplay, 'distance');
% compute pixel size in visual degrees
displayPixelSizeDegrees = ...
    2 * atand(0.5*displayPixelSizeMeters/viewingDistanceMeters);
% divide the stimulus size in degrees by the pixel size in degrees to get 
% the number of stimulus pixels
stimParams.pixelsAlongWidthDim = ...
    round(stimParams.sizeDegs/displayPixelSizeDegrees(1));
stimParams.pixelsAlongHeightDim = ...
    round(stimParams.sizeDegs/displayPixelSizeDegrees(2));

% generate the scene

scene = generateGaborScene('stimParams', stimParams);
% Place the scene at the same distance as the display
scene = sceneSet(scene, 'distance', viewingDistanceMeters);
% Give the scene a name
scene = sceneSet(scene, 'name', 'GABOR SCENE');

%% Visualizing the scene
visualizeScene(scene);

%% Creating the scene on display

% Convert the idealized stimulus scene into a scene that represents 
% the stimulus as realized by the presentation display
displayedScene = realizeSceneOnDisplay(scene, presentationDisplay);
% Give the scene a name
displayedScene = sceneSet(displayedScene, 'name', 'GABOR SCENE (displayed)');

%% Visualizing the display scene
visualizeScene(displayedScene);










