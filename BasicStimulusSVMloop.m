accuracy_FEM = [];
accuracy_noFEM = [];
for l = 0:1:12
    % Stimulus 1 
    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', 10, ... 
        'orientationDegs', 90, ...             
        'phaseDegs', 180, ...                    
        'sizeDegs', 0.5, ...                    
        'sigmaDegs', 0.2/3, ...                 
        'contrast', l/100,...                  
        'meanLuminanceCdPerM2', 15, ...         
        'pixelsAlongWidthDim', [], ...          
        'pixelsAlongHeightDim', [] ...          
        );
    scene1 = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay, ...
        'minimumPixelsPerHalfPeriod', 10);
    % Stimulus 2
    stimParams.contrast = 0;
    scene2 = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay, ...
        'minimumPixelsPerHalfPeriod', 10);
    % visualizeScene(scene1, ...
    %     'displayRadianceMaps', false);
    % visualizeScene(scene2, ...
    %     'displayRadianceMaps', false);

    % Generate wavefront-aberration derived human optics
    theOI = oiCreate('wvf human');
    OI_1 = oiCompute(theOI, scene1);
    OI_2 = oiCompute(theOI, scene2);
    % visualizeOpticalImage(OI_1, ...
    %     'displayRadianceMaps', false);
    % visualizeOpticalImage(OI_2, ...
    %     'displayRadianceMaps', false);

    integrationTime = 50;
    % Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
    theMosaic = cMosaic('sizeDegs', [1, 1] * stimParams.sizeDegs, ... 
                'integrationTime', integrationTime / 1000); % 10 msec integration time
    % with eye movement
    duration = 500 / 1000;
    nTrials = (duration*1000) / integrationTime;
    nTrialsNum = 50;
    noisyExcitationResponseInstances1 = zeros(nTrialsNum, nTrials, 3822);
    noisyExcitationResponseInstances2 = zeros(nTrialsNum, nTrials, 3822);
    for i = 1:nTrialsNum
        theMosaic.emGenSequence(duration, 'nTrials', nTrials, 'centerPaths', true,'microsaccadeType', 'none');
        [~,noisyExcitationResponseInstances_11, ~,~,~] = theMosaic.compute(OI_1, 'nTrials', 1, 'withFixationalEyeMovements', true);
        noisyExcitationResponseInstances1(i,:,:) = noisyExcitationResponseInstances_11;
        [~,noisyExcitationResponseInstances_22, ~,~,~] = theMosaic.compute(OI_2, 'nTrials', 1, 'withFixationalEyeMovements', true);
        noisyExcitationResponseInstances2(i,:,:) = noisyExcitationResponseInstances_22;
    end

    % without eye movement
    [noiseFreeExcitationResponseInstances1_noFEM,noisyExcitationResponseInstances1_noFEM, ~,~,~] = theMosaic.compute(OI_1, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
    [noiseFreeExcitationResponseInstances2_noFEM,noisyExcitationResponseInstances2_noFEM, ~,~,~] = theMosaic.compute(OI_2, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);

    % Simulate a 2AFC task 
    taskIntervals = 2;
    [classificationMatrix, classLabels] = generateSetUpForClassifier(...
        noisyExcitationResponseInstances1, noisyExcitationResponseInstances2, taskIntervals);

    % Find principal components of the responses
    [pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

    % Project the responses onto the space formed by the first 4 PC vectors
    pcComponentsNumForClassification = 2;
    classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);

    % Visualize the classification matrix and its projection to the PC space
    % visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
    svm = fitcsvm(classificationMatrixProjection,classLabels);
    % Measure performance of the SVM classifier using a 10-fold crossvalidation
    % approach
    % Perform a 10-fold cross-validation on the trained SVM model
    kFold = 10;
    CVSVM = crossval(svm,'KFold',kFold);

    % Compute classification loss for the in-sample responses using a model 
    % trained on out-of-sample responses
    fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
    % Average percent correct across all folds 
    percentCorrect_FEM = mean(fractionCorrect)*100;
    accuracy_FEM(end+1) = percentCorrect_FEM;

    % Simulate a 2AFC task 
    taskIntervals = 2;
    [classificationMatrix, classLabels] = generateSetUpForClassifier_noEM(...
        noisyExcitationResponseInstances1_noFEM, noisyExcitationResponseInstances2_noFEM, taskIntervals);

    % Find principal components of the responses
    [pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

    % Project the responses onto the space formed by the first 4 PC vectors
    pcComponentsNumForClassification = 2;
    classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);

    % Visualize the classification matrix and its projection to the PC space
    % visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
    % Train a binary SVM classifier and visualize the support vectors in 2
    % dimensions
    svm = fitcsvm(classificationMatrixProjection,classLabels);
    % Visualize the data along with the hyperplane computed by the SVM 
    % visualizeSVMmodel(svm, classificationMatrixProjection, classLabels);
    % Measure performance of the SVM classifier using a 10-fold crossvalidation approach
    % Perform a 10-fold cross-validation on the trained SVM model
    kFold = 10;
    CVSVM = crossval(svm,'KFold',kFold);

    % Compute classification loss for the in-sample responses using a model 
    % trained on out-of-sample responses
    fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
    % Average percent correct across all folds 
    percentCorrect_noFEM = mean(fractionCorrect)*100;
    accuracy_noFEM(end+1) = percentCorrect_noFEM;
end
%% figure
figure
x = [0:1:12];
y1 = accuracy_noFEM;
y2 = accuracy_FEM;
y1_up = y1 + std(y1);
y1_down = y1 - std(y1);
patch([x fliplr(x)], [y1_up  fliplr(y1_down)], 'r','FaceAlpha',0.2, 'EdgeColor','none'); 
hold on;
plot(x, y1, '-or','LineWidth', 2);
y2_up = y2 + std(y2);
y2_down = y2 - std(y2);
ylim([0 100])
patch([x fliplr(x)], [y2_up  fliplr(y2_down)], 'b','FaceAlpha',0.2, 'EdgeColor','none'); 
plot(x, y2, '-ob', 'LineWidth', 2);
title('Changing in SVM Acuracy (luminance = 15)','fontsize', 15)
xlabel('Stimulus Contrast (%)')
ylabel('Accuracy (%)')
legend('','no eye movement','','with eye movement')
hold off

%% Supporting function
function [classificationMatrix, classLabels] = generateSetUpForClassifier(coneExcitationsTest, coneExcitationsNull, taskIntervals)
    nTrials = size(coneExcitationsTest, 1);
    responseSize = size(coneExcitationsTest, 3)*size(coneExcitationsTest, 2);
    testResponses = reshape(coneExcitationsTest,nTrials,responseSize);
    nullResponses = reshape(coneExcitationsNull,nTrials,responseSize);
        
    % Assemble the response vectors into a classification matrix simulating either 
    % a one interval task or a two-interval task.
    if (taskIntervals == 1)
        % In the one interval task, the null and test response instances are labelled as the 2 classes.
        % Allocate matrices
        classificationMatrix = nan(2*nTrials, responseSize);
        classLabels = nan(2*nTrials, 1);
        % Class 1
        classificationMatrix(1:nTrials,:) = nullResponses;
        classLabels((1:nTrials)) = 0;
        % Class 2
        classificationMatrix(nTrials+(1:nTrials),:) = testResponses;
        classLabels(nTrials+(1:nTrials)) = 1;
    elseif (taskIntervals == 2)
        % In the two inteval task, we concatenate [null test] as one class and [test null] as the other. 
        % Allocate matrices
        classificationMatrix = nan(nTrials, 2*responseSize);
        classLabels = nan(nTrials, 1);
        halfTrials = floor(nTrials/2);
        % Class 1
        classificationMatrix(1:halfTrials,:) = [...
            nullResponses(1:halfTrials,:) ...
            testResponses(1:halfTrials,:)];
        classLabels((1:halfTrials)) = 0;
        % Class 2
        idx = halfTrials+(1:halfTrials);
        classificationMatrix(idx,:) = [...
            testResponses(idx,:) ...
            nullResponses(idx,:)];
        classLabels(idx) = 1;
    else
        error('Task can have 1 or 2 intervals only.')
    end
end

function [classificationMatrix, classLabels] = generateSetUpForClassifier_noEM(coneExcitationsTest, coneExcitationsNull, taskIntervals)
    nTrials = size(coneExcitationsTest, 1);
    responseSize = size(coneExcitationsTest, 3);
    testResponses = squeeze(coneExcitationsTest);
    nullResponses = squeeze(coneExcitationsNull);
        
    % Assemble the response vectors into a classification matrix simulating either 
    % a one interval task or a two-interval task.
    if (taskIntervals == 1)
        % In the one interval task, the null and test response instances are labelled as the 2 classes.
        % Allocate matrices
        classificationMatrix = nan(2*nTrials, responseSize);
        classLabels = nan(2*nTrials, 1);
        % Class 1
        classificationMatrix(1:nTrials,:) = nullResponses;
        classLabels((1:nTrials)) = 0;
        % Class 2
        classificationMatrix(nTrials+(1:nTrials),:) = testResponses;
        classLabels(nTrials+(1:nTrials)) = 1;
    elseif (taskIntervals == 2)
        % In the two inteval task, we concatenate [null test] as one class and [test null] as the other. 
        % Allocate matrices
        classificationMatrix = nan(nTrials, 2*responseSize);
        classLabels = nan(nTrials, 1);
        halfTrials = floor(nTrials/2);
        % Class 1
        classificationMatrix(1:halfTrials,:) = [...
            nullResponses(1:halfTrials,:) ...
            testResponses(1:halfTrials,:)];
        classLabels((1:halfTrials)) = 0;
        % Class 2
        idx = halfTrials+(1:halfTrials);
        classificationMatrix(idx,:) = [...
            testResponses(idx,:) ...
            nullResponses(idx,:)];
        classLabels(idx) = 1;
    else
        error('Task can have 1 or 2 intervals only.')
    end
end


function visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
    figure(); clf;
    subplot(2,1,1);
    visualizeClassificationMatrix(classificationMatrix, taskIntervals, 'cone index');
    subplot(2,1,2);
    visualizeClassificationMatrix(classificationMatrixProjection, [], 'principal component no');
end

function visualizeClassificationMatrix(classificationMatrix, taskIntervals, xAxisLabel)
    imagesc(1:size(classificationMatrix,2), 1:size(classificationMatrix,1), classificationMatrix);
    if (~isempty(taskIntervals))
        hold on
        if taskIntervals == 1
            %plot(size(classificationMatrix,2)/2*[1 1], [1 size(classificationMatrix,1)], 'c-', 'LineWidth', 2.0);
            plot([1 size(classificationMatrix,2)], size(classificationMatrix,1)/2*[1 1], 'c-', 'LineWidth', 2.0);
        else
            plot(size(classificationMatrix,2)/2*[1 1], [1 size(classificationMatrix,1)], 'c-', 'LineWidth', 2.0);
            set(gca, 'XTick', [])
        end
        hold off
    end
    axis 'xy'
    set(gca, 'FontSize',14)
    xlabel(xAxisLabel)
    ylabel('trials')
    if (strcmp(xAxisLabel, 'principal component no'))
        set(gca, 'XTick', 1:size(classificationMatrix,2))
    end
    title('classification matrix')
    colormap(gray)
end

function visualizePrincipalComponents(pcVectors, varianceExplained, theMosaic)
    figure(); clf;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 2, ...
       'heightMargin',  0.1, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.01, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.1, ...
       'topMargin',      0.1);
    for pcaComponentIndex = 1:4
        title = sprintf('PCA %d, variance\nexplained: %2.2f%%', ...
            pcaComponentIndex, varianceExplained(pcaComponentIndex));
        r = floor((pcaComponentIndex-1)/2)+1;
        c = mod(pcaComponentIndex-1,2)+1;
        ax = subplot('Position', subplotPosVectors(r,c).v);
        theMosaic.visualize('activation', pcVectors(:,pcaComponentIndex))        
        ylabel('')
        set(ax, 'YTick', [])
        if (pcaComponentIndex < 3)
            xlabel('')
            set(ax, 'XTick', [])
        end
    end
end

function visualizeSVMmodel(svmModel, data, classes)
    sv = svmModel.SupportVectors;
    h = max(abs(data(:)))/100; % Mesh grid step size
    r = -h*100:h:h*100;
    [X1,X2] = ndgrid(r, r);
    [~,score] = predict(svmModel,[X1(:),X2(:)]);
    scoreGrid = reshape(score(:,1),numel(r), numel(r));

    figure(); clf;
    class0Indices = find(classes == 0);
    class1Indices = find(classes == 1);
    
    plot(data(class0Indices,1),data(class0Indices,2),'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'c'); hold on
    plot(data(class1Indices,1), data(class1Indices,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    contourf(X1,X2,scoreGrid, 50); 
    plot(data(class0Indices,1),data(class0Indices,2),'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'c')
    plot(data(class1Indices,1), data(class1Indices,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(sv(:,1),sv(:,2),'ks','MarkerSize',12, 'LineWidth', 1.5);
    colormap(brewermap(1024, 'RdBu'))
    hold off
    xlabel('PC component #1 activation')
    ylabel('PC component #2 activation')
    legend('stimulus 1', 'stimulus 2')
    set(gca, 'FontSize',14)
    axis 'square'
end

