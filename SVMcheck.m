stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 10, ... 
    'orientationDegs', 90, ...               
    'phaseDegs', 90, ...                    
    'sizeDegs', 0.5, ...                    
    'sigmaDegs', 0.25/3, ...                 
    'contrast', 0/100,...                  
    'meanLuminanceCdPerM2', 15, ...
    'center', [0 0], ...
    'pixelsAlongWidthDim', [], ...          
    'pixelsAlongHeightDim', [] ...          
    );
% get response by the function computeConeResponseforSVM
ResponsesTest = computeConeResponseforSVM(stimParams, 'null', 'contrast', 0, 'nTrials', 10, 'responseFlag', 'excitation');
ResponsesNull = computeConeResponseforSVM(stimParams, 'null', 'contrast', 0, 'nTrials', 10, 'responseFlag', 'excitation');
noisyResponsesTest = ResponsesTest.noisyExcitation;
noisyResponsesNull = ResponsesNull.noisyExcitation;
% Simulate a 2AFC task 
taskIntervals = 2;
[classificationMatrix, classLabels] = generateSetUpForClassifier(...
    noisyResponsesTest(:,1,:), noisyResponsesNull(:,1,:), taskIntervals, 'false');

% Find principal components of the responses
[pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

% Project the responses onto the space formed by the first 4 PC vectors
pcComponentsNumForClassification = 2;
classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);

% Visualize the classification matrix and its projection to the PC space
visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
% Train a binary SVM classifier and visualize the support vectors in 2 dimensions
svm = fitcsvm(classificationMatrixProjection,classLabels);

% Visualize the data along with the hyperplane computed by the SVM 
visualizeSVMmodel(svm, classificationMatrixProjection, classLabels);

% Measure performance of the SVM classifier using a 10-fold crossvalidation approach
% Perform a 10-fold cross-validation on the trained SVM model
CVSVM = crossval(svm,'KFold',10);

% Compute classification loss for the in-sample responses using a model 
% trained on out-of-sample responses
fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
% Average percent correct across all folds 
percentCorrect = mean(fractionCorrect)*100
noisyAccuracy(:,count) = fractionCorrect*100;

disp(percentCorrect)