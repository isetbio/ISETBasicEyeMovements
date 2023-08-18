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
    t = tiledlayout('flow');
    t.TileSpacing = 'tight';
    title(t,"Responses of cone index " + targetConeID,'FontWeight','bold')
    % xlabel(t,'time (seconds)','FontSize',14)
    
    nexttile
    plot(timeAxis, squeeze(noisyResponses_ex(:,:,targetConeID)), 'c-', 'color', '#7F7FD3');
    hold on;
    plot(timeAxis, squeeze(mean(noisyResponses_ex(:,:,targetConeID),1)), 'LineWidth', 2.5, 'DisplayName','Excitation (R*/cone/tau)', 'color', '#0000A7');
    plot(timeAxis, squeeze(mean(noiseFreeResponses_ex(:,:,targetConeID),1)), 'k-', 'LineWidth', 1.5);
    ylabel('Excitation (R*/cone/tau)');
    set(gca,'linewidth',1)
    set(gca,'FontSize', 12)
    hold off;
    if ~isempty(excitationScale)
        ylim(excitationScale)
    end
    % title("Responses of cone index " + targetConeID)
    
    nexttile
    plot(timeAxis, squeeze(noisyResponses_ph(:,:,targetConeID)), 'm-', 'color', '#D37FD3');
    hold on;
    plot(timeAxis, squeeze(mean(noisyResponses_ph(:,:,targetConeID),1)), 'r-', 'LineWidth', 2.5, 'color', '#A700A7');
    % plot(timeAxis, squeeze(mean(noisyResponses_ph(:,:,targetConeID),1)), 'r-', 'LineWidth', 2.5, 'DisplayName','Photocurrent (pAmps)');
    plot(timeAxis, squeeze(mean(noiseFreeResponses_ph(:,:,targetConeID),1)), 'k-', 'LineWidth', 1.5);
    hold off;
    ylabel('Photocurrent (pAmps)');
    xlabel('time (seconds)','FontSize',15);
    set(gca,'linewidth',1)
    set(gca,'FontSize', 12)
    if ~isempty(photocurrentScale)
        ylim(photocurrentScale)
    end
    
    meanEx = mean(noiseFreeResponses_ex(:,:,targetConeID));
    meanPhoto = mean(noiseFreeResponses_ph(:,:,targetConeID));
end