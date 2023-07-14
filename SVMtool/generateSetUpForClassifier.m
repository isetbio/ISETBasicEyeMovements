function [classificationMatrix, classLabels] = generateSetUpForClassifier(coneResponseTest, coneResponseNull, taskIntervals, hasFEM)
    if notDefined('hasFEM'), hasFEM = 'true'; end
    nTrials = size(coneResponseTest, 1);
    switch hasFEM
        case 'true'
            responseSize = size(coneResponseTest, 3)*size(coneResponseTest, 2);
            testResponses = reshape(coneResponseTest,nTrials,responseSize);
            nullResponses = reshape(coneResponseNull,nTrials,responseSize);
        case 'false'
            responseSize = size(coneResponseTest, 3);
            testResponses = squeeze(coneResponseTest);
            nullResponses = squeeze(coneResponseNull);
        otherwise
            error('Unknown FEM info.');
    end
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
    