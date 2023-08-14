function [meanEx, meanPhoto] = visualizeAllResponses(responseStruct, varargin)
    p = inputParser;
    p.addParameter('targetCone', '', @(x) (isempty(x) | isnumeric(x)));
    p.addParameter('excitationScale', '');
    p.addParameter('photocurrentScale', '');
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    targetCone = p.Results.targetCone;
    excitationScale = p.Results.excitationScale;
    photocurrentScale = p.Results.photocurrentScale;
    
    noiseFreeResponses_ex = responseStruct.noiseFreeExcitation;
    noisyResponses_ex = responseStruct.noisyExcitation;
    noiseFreeResponses_ph = responseStruct.noiseFreePhotocurr;
    noisyResponses_ph = responseStruct.noisyPhotocurr;
    timeAxis = responseStruct.timeAxis;
    
    [~,idx] = max(noiseFreeResponses_ex(:));
    [~,~,maxConeID] = ind2sub(size(noiseFreeResponses_ex), idx);
    if isempty(targetCone)
        targetConeID = maxConeID;
    else
        targetConeID = targetCone;
    end
    
    figure()
    % set(gca, 'FontSize', 10);
    tiledlayout('flow')
    nexttile
    plot(timeAxis, squeeze(noisyResponses_ph(:,:,targetConeID)), 'm-', 'DisplayName','');
    hold on;
    plot(timeAxis, squeeze(mean(noisyResponses_ph(:,:,targetConeID),1)), 'r-', 'LineWidth', 2.5, 'DisplayName','Photocurrent (pAmps)');
    plot(timeAxis, squeeze(mean(noiseFreeResponses_ph(:,:,targetConeID),1)), 'k-', 'LineWidth', 1.5);
    hold off;
    ylabel('Photocurrent (pAmps)');
    if ~isempty(photocurrentScale)
        ylim(photocurrentScale)
    end
    title("Responses of cone index " + targetConeID)
    
    nexttile
    plot(timeAxis, squeeze(noisyResponses_ex(:,:,targetConeID)), 'c-', 'DisplayName','');
    hold on;
    plot(timeAxis, squeeze(mean(noisyResponses_ex(:,:,targetConeID),1)), 'b-', 'LineWidth', 2.5, 'DisplayName','Excitation (R*/cone/tau)');
    plot(timeAxis, squeeze(mean(noiseFreeResponses_ex(:,:,targetConeID),1)), 'k-', 'LineWidth', 1.5);
    ylabel('Excitation (R*/cone/tau)');
    hold off;
    xlabel('time (seconds)');
    if ~isempty(excitationScale)
        ylim(excitationScale)
    end
    
    meanEx = mean(noiseFreeResponses_ex(:,:,targetConeID));
    meanPhoto = mean(noiseFreeResponses_ph(:,:,targetConeID));
end