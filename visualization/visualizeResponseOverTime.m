function [meanNoiseFree, meanNoisy, coneID] = visualizeResponseOverTime(responseStruct, responseType, varargin)
    p = inputParser;
    p.addParameter('targetCone', @(x) (isempty(x) | isnumeric(x)));
%     p.addParameter('auto', 'true', @(x) (isempty(x) | ischar(x)));
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    targetCone = p.Results.targetCone;
%     auto = p.Results.auto;
    
    figure()
    if strcmp(responseType,'excitation')
        noiseFreeResponses = responseStruct.noiseFreeExcitation;
        noisyResponses = responseStruct.noisyExcitation;
        ylabel('excitations per integration time (R*/cone/tau)');
    elseif strcmp(responseType,'photocurrent')
        noiseFreeResponses = responseStruct.noiseFreePhotocurr;
        noisyResponses = responseStruct.noisyPhotocurr;
        ylabel('Photocurrent (pAmps)');
    end
    [~,idx] = max(noiseFreeResponses(:));
    [~,~,maxConeID] = ind2sub(size(noiseFreeResponses), idx);
    if isempty(targetCone)
        targetConeID = maxConeID;
    else
        targetConeID = targetCone;
    end
    % Plot the time series response for individual instances
    timeAxis = responseStruct.timeAxis;
    plot(timeAxis, squeeze(noisyResponses(:,:,targetConeID)), '.', 'color', [.5 .5 .5]);
    hold on;
    % Plot the time series response for the mean of the individual instances
    plot(timeAxis, squeeze(mean(noisyResponses(:,:,targetConeID),1)), 'g-', 'LineWidth', 2.0);
    % Plot the noise-free time series response in red
    plot(timeAxis, squeeze(mean(noiseFreeResponses(:,:,targetConeID),1)), 'r', 'LineWidth', 1.5);
    xlabel('time (seconds)');
    set(gca, 'FontSize', 16);
    meanNoiseFree = mean(noiseFreeResponses(:,:,targetConeID));
    meanNoisy = mean(noisyResponses(:,:,targetConeID));
    coneID = targetConeID
end