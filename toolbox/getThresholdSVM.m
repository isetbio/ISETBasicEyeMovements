function [threshold, curvefit, gof, data, original_data] = getThresholdSVM(stimParamStruct, theMosaic, responseType, varargin)

% Compute function for detection threshold of different grating contrast
%
% Syntax:
%   [threshold, curvefit, gof] = getThresholdSVM(stimParamStruct,responseType, varargin); 
%   if ignoring curvefit and gof, just call
%   [threshold, ~ ~] = getThresholdSVM(stimParamStruct, responseType, varargin)
%
% Inputs:
%    stimParamStruct                - a struct containing properties of the
%                                     the grating.
%    responseType                   - Specify expected response type. valid string are
%                                       - 'excitation'
%                                       - 'photocurrent'
%
% Optional key/value input arguments:
%    'visualize'                   - whether to visualize the fitted psychometric curve. 
%                                    Default: false
%                                    Valid values are: 
%                                        - 'true' (visualize the curve)
%                                        - 'false' 
%    'showResult'                   - whether to print out the threshold after calculation. 
%                                     Default: false
%                                     Valid values are: 
%                                        - 'true' (print out on command window)
%                                        - 'false' 
%    'nTrials'                      - Integer. Specify the number of trials
%                                     created for each time point
%    'center'                       - [x y] of the center for the grating. 
%                                     default [0 0]
%    'kFold'                        - Integer. Specify the number of fold
%                                     to split for cross-validation 
%                                     defualt 10
%    'min'                          - Integer. minimum number to start the
%                                     looping. Default: 0
%    'nstep'                        - Integer. Specifies how big the each step of the loop
%                                     Default: 1
%    'max'                          - Integer. maximum number to start the
%                                     looping. Default: 10
%
% Outputs:
%    threshold  - A value in percentage specifies the detection threshold of a given stimulus. 
%
%    curvefit   - A cfit object that specify its model used to fit the
%                 curve and the coefficient of its auto-bestfit model
%
%    gof        - A Struct that contain info about the Goodness of Fit of the fitted curve. 
%                 access specific gof by calling:
%                   - gof.sse
%                   - gof.rsquare
%                   - gof.dfe
%                   - gof.adjrsquare
%                   - gof.rmse
%
% History:
%    13/07/23  Qiyuan Feng

    p = inputParser;
    p.addParameter('visualize', '', @(x) (isempty(x) | ischar(x)));
    p.addParameter('kFold', 10, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('min', 0, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('nstep', 1, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('max', 10, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('nTrials', 10, @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('center', [0 0]);
    p.addParameter('showResult', '', @(x) (isempty(x) | ischar(x)));
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    visualize = p.Results.visualize;
    kFold = p.Results.kFold;
    min = p.Results.min;
    nstep = p.Results.nstep;
    max = p.Results.max;
    nTrials = p.Results.nTrials;
    center = p.Results.center;
    showResult = p.Results.showResult;
    
    accuracy = [];
    st_deviation = [];
    original = zeros(kFold, ((max-min)/nstep));
    count = 1;

    for c = min:nstep:max
        % get response by the function computeConeResponseforSVM
        if strcmp(responseType,'excitation')
            ResponsesTest = computeConeResponseforSVM(stimParamStruct, theMosaic, 'test', 'contrast', c, 'nTrials', nTrials, 'center', center, 'responseFlag', 'excitation');
            ResponsesNull = computeConeResponseforSVM(stimParamStruct, theMosaic, 'null', 'contrast', 0, 'nTrials', nTrials, 'responseFlag', 'excitation');
            noisyResponsesTest = ResponsesTest.noisyExcitation;
            noisyResponsesNull = ResponsesNull.noisyExcitation;
        elseif strcmp(responseType,'photocurrent')
            ResponsesTest = computeConeResponseforSVM(stimParamStruct, theMosaic, 'test', 'contrast', c, 'nTrials', nTrials, 'center', center);
            ResponsesNull = computeConeResponseforSVM(stimParamStruct, theMosaic, 'null', 'nTrials', nTrials);
            noisyResponsesTest = ResponsesTest.noisyPhotocurr;
            noisyResponsesNull = ResponsesNull.noisyPhotocurr;
        end
        % Simulate a 2AFC task 
        taskIntervals = 2;
        [classificationMatrix, classLabels] = generateSetUpForClassifier(...
            noisyResponsesTest, noisyResponsesNull, taskIntervals, 'true');

        % Find principal components of the responses
        [pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

        % Project the responses onto the space formed by the first 4 PC vectors
        pcComponentsNumForClassification = 2;
        classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);

        if strcmp(visualize, 'matrix') || strcmp(visualize, 'all')
            % Visualize the classification matrix and its projection to the PC space
            visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
        end
        % Train a binary SVM classifier and visualize the support vectors in 2 dimensions
        svm = fitcsvm(classificationMatrixProjection,classLabels);

        if strcmp(visualize, 'SVM') || strcmp(visualize, 'all')
            % Visualize the data along with the hyperplane computed by the SVM 
            visualizeSVMmodel(svm, classificationMatrixProjection, classLabels);
        end

        % Measure performance of the SVM classifier using a 10-fold crossvalidation approach
        % Perform a 10-fold cross-validation on the trained SVM model
        CVSVM = crossval(svm,'KFold',kFold);

        % Compute classification loss for the in-sample responses using a model 
        % trained on out-of-sample responses
        fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
        % Average percent correct across all folds 
        percentCorrect = mean(fractionCorrect)*100;
        accuracy(end+1) = percentCorrect;
        original(:,count) = fractionCorrect*100;
        st_deviation(end+1) = std(original(:,count));
        count = count +1;
    end
    
    
    x = [min:nstep:max];
    [curvefit, gof] = fitLogistic(x, accuracy); % validate the rsquare, making sure it's > 0.9
    sse = gof.sse;
    rsquare = gof.rsquare;
    setThreshold = 75;
    find = @(x) setThreshold - curvefit(x);
    threshold = fzero(find, 0);
    
    if strcmp(showResult,'true')
        fprintf('Detection Threshold is %.4f with a r-square of %.4f\n', threshold, rsquare);
    end
    
    if strcmp(visualize,'true')
        figure()
        plot(x,accuracy,'.', 'MarkerSize', 12);
        hold on;
        bestfit = plot(curvefit);
        set(bestfit, 'LineWidth', 2.5, 'color', [0.4940 0.1840 0.5560])
        plot([0 threshold], [setThreshold setThreshold], '--','color', [0.5 0.5 0.5],'LineWidth', 2);
        plot([threshold threshold], [0 setThreshold], '--','color', [0.5 0.5 0.5],'LineWidth', 2);
        legend('original data','fitted curve')
        xlabel('Stimulus Contrast (%)','fontsize', 12)
        ylabel('Accuracy (%)','fontsize', 12)
        ylim([0 101.5])
        if ~isnan(threshold)
            title({
            ["Threshold = " + threshold + "%"]
            ["SSE: " + sse + ", R-square: " + rsquare]
            }, 'fontsize', 15);
        end
        hold off;
    end
    data = accuracy;
    original_data = original;
end

function [curvefit, gof] = fitLogistic(x,y)
    [curvefit, gof] = fit(x(:), y(:), 'a+b/(1+exp(-c*(x-d)))');
    iteration = 1;
    while gof.rsquare < 0.9 && iteration < 5
        fprintf('poor curve fitting, trying again.\n');
        [curvefit, gof] = fit(x(:), y(:), 'a+b/(1+exp(-c*(x-d)))');
        iteration = iteration + 1;
    end
    if gof.rsquare < 0.9
        fprintf('Warning: poor auto curve fitting after 5 iterations\n');
    end
end
