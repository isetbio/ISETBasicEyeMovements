accuracy_FEM = [];
accuracy_noFEM = [];
std_noFEM = [];
kFold = 10;
nmin = 0;
nstep = 0.5;
nmax = 15;
noisyAccuracy_noFEM = zeros(kFold, ((nmax-nmin)/nstep));
count = 1;
nTrialsNum = 100;
for l = nmin:nstep:nmax
    % Stimulus 1 
    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', 10, ... 
        'orientationDegs', 90, ...             
        'phaseDegs', 0, ...                    
        'sizeDegs', 0.5, ...                    
        'sigmaDegs', 0.25/3, ...                 
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
%     % with eye movement
%     duration = 500 / 1000;
%     nTrials = (duration*1000) / integrationTime;
%     nTrialsNum = 50;
%     noisyExcitationResponseInstances1 = zeros(nTrialsNum, nTrials, 3822);
%     noisyExcitationResponseInstances2 = zeros(nTrialsNum, nTrials, 3822);
%     for i = 1:nTrialsNum
%         theMosaic.emGenSequence(duration, 'nTrials', nTrials, 'centerPaths', true,'microsaccadeType', 'none');
%         [~,noisyExcitationResponseInstances_11, ~,~,~] = theMosaic.compute(OI_1, 'nTrials', 1, 'withFixationalEyeMovements', true);
%         noisyExcitationResponseInstances1(i,:,:) = noisyExcitationResponseInstances_11;
%         [~,noisyExcitationResponseInstances_22, ~,~,~] = theMosaic.compute(OI_2, 'nTrials', 1, 'withFixationalEyeMovements', true);
%         noisyExcitationResponseInstances2(i,:,:) = noisyExcitationResponseInstances_22;
%     end

    % without eye movement
    [noiseFreeExcitationResponseInstances1_noFEM,noisyExcitationResponseInstances1_noFEM, ~,~,~] = theMosaic.compute(OI_1, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);
    [noiseFreeExcitationResponseInstances2_noFEM,noisyExcitationResponseInstances2_noFEM, ~,~,~] = theMosaic.compute(OI_2, 'nTrials', nTrialsNum, 'withFixationalEyeMovements', false);

%     % Simulate a 2AFC task 
%     taskIntervals = 2;
%     [classificationMatrix, classLabels] = generateSetUpForClassifier(...
%         noisyExcitationResponseInstances1, noisyExcitationResponseInstances2, taskIntervals);
% 
%     % Find principal components of the responses
%     [pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);
% 
%     % Project the responses onto the space formed by the first 4 PC vectors
%     pcComponentsNumForClassification = 2;
%     classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);
% 
%     % Visualize the classification matrix and its projection to the PC space
%     % visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
%     svm = fitcsvm(classificationMatrixProjection,classLabels);
%     % Measure performance of the SVM classifier using a 10-fold crossvalidation
%     % approach
%     % Perform a 10-fold cross-validation on the trained SVM model
%     kFold = 10;
%     CVSVM = crossval(svm,'KFold',kFold);
% 
%     % Compute classification loss for the in-sample responses using a model 
%     % trained on out-of-sample responses
%     fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
%     % Average percent correct across all folds 
%     percentCorrect_FEM = mean(fractionCorrect)*100;
%     accuracy_FEM(end+1) = percentCorrect_FEM;

    % Simulate a 2AFC task 
    taskIntervals = 2;
    [classificationMatrix, classLabels] = generateSetUpForClassifier(...
        noisyExcitationResponseInstances1_noFEM, noisyExcitationResponseInstances2_noFEM, taskIntervals, 'false');

    % Find principal components of the responses
    [pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

    % Project the responses onto the space formed by the first 4 PC vectors
    pcComponentsNumForClassification = 2;
    classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);

    % Visualize the classification matrix and its projection to the PC space
    % visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
    % Train a binary SVM classifier and visualize the support vectors in 2 dimensions
    svm = fitcsvm(classificationMatrixProjection,classLabels);
    % Visualize the data along with the hyperplane computed by the SVM 
    % visualizeSVMmodel(svm, classificationMatrixProjection, classLabels);
    % Measure performance of the SVM classifier using a 10-fold crossvalidation approach
    % Perform a 10-fold cross-validation on the trained SVM model
    CVSVM = crossval(svm,'KFold',kFold);

    % Compute classification loss for the in-sample responses using a model 
    % trained on out-of-sample responses
    fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
    % Average percent correct across all folds 
    percentCorrect_noFEM = mean(fractionCorrect)*100;
    accuracy_noFEM(end+1) = percentCorrect_noFEM;
    noisyAccuracy_noFEM(:,count) = fractionCorrect*100;
    std_noFEM(end+1) = std(noisyAccuracy_noFEM(:,count));
    count = count +1;
end
%% figure
figure
x = [nmin:nstep:nmax];
y1 = accuracy_noFEM;
% y2 = accuracy_FEM;
y1_up = y1 + std_noFEM;
y1_down = y1 - std_noFEM;
% errorbar(x,y1,std_noFEM, '-ob', 'LineWidth', 2)
patch([x fliplr(x)], [y1_up  fliplr(y1_down)], 'r','FaceAlpha',0.2, 'EdgeColor','none'); 
hold on;
plot(x, y1, '-or','LineWidth', 2);
% y2_up = y2 + std(y2);
% y2_down = y2 - std(y2);
% ylim([0 100])
% patch([x fliplr(x)], [y2_up  fliplr(y2_down)], 'b','FaceAlpha',0.2, 'EdgeColor','none'); 
% plot(x, y2, '-ob', 'LineWidth', 2);
title('Changing in SVM Acuracy (luminance = 15)','fontsize', 15)
xticks(nmin:1:nmax);
ylim([0 100])
yticks(0:10:100);
xlabel('Stimulus Contrast (%)')
ylabel('Accuracy (%)')
% legend('','no eye movement','','with eye movement')
hold off


% boxplot(noisyAccuracy_noFEM)
% title('Changing in SVM Acuracy (luminance = 15)','fontsize', 15)
% xlabel('Stimulus Contrast (%)')
% ylabel('Accuracy (%)')
% legend('','no eye movement','','with eye movement')

%% Try curve fitting
% paramfit = wblfit(accuracy_noFEM);
% scale = paramfit(1);
% shape = paramfit(2);
% plot(x,wblcdf(accuracy_noFEM,scale,shape),'DisplayName','A=10, B=2')

plot(x,accuracy_noFEM,'o');
hold on;
% modelFun =  @(p,x) p(3) .* (x./p(1)).^(p(2)-1) .* exp(-(x./p(1)).^p(2));
% startingVals = [10 2 5];
% nlModel = fitnlm(x,accuracy_noFEM,modelFun,startingVals);
% line(x,predict(nlModel,x),'Color','r');

[curvefit, gof] = fit(x(:), accuracy_noFEM(:), 'a+b/(1+exp(-c*(x-d)))');
% [curvefit, gof] = fit(x(:), accuracy_noFEM(:), 'c*a*b*x^(b-1)*exp(-a*x^b)', 'StartPoint', [0.01, 2, 5]);
sse = gof.sse;
rsquare = gof.rsquare
plot(curvefit)
hold off;


