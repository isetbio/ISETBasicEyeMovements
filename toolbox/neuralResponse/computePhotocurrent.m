function photocurrentResponses = computePhotocurrent(coneExcitationResponses, temporalSupportSeconds, noiseFlag, cutoffFlag)
    if notDefined('cutoffFlag'), cutoffFlag = 'false'; end
    % Compute photocurrent responses from the cone excitation responses
    photocurrentResponses = 0*coneExcitationResponses;
    switch cutoffFlag
        case 'true'
            impulseResponse = computePhotocurrentImpulseResponse(temporalSupportSeconds, 'true');
        case 'false'
            impulseResponse = computePhotocurrentImpulseResponse(temporalSupportSeconds, 'false');
    end
    impulseResponse(1) = 0.0;
    dt = temporalSupportSeconds(2)-temporalSupportSeconds(1);
    
    for iTrial = 1:size(coneExcitationResponses,1)
        for iCone = 1:size(coneExcitationResponses,3)
            coneExcitation = squeeze(coneExcitationResponses(iTrial,:,iCone));
            pCurrent = conv(coneExcitation, impulseResponse, 'full');
            pCurrent = pCurrent(1:length(coneExcitation));
            if (~strcmp(noiseFlag, 'none'))
                photocurrentResponses(iTrial,:,iCone) = osAddNoise(pCurrent, 'sampTime', dt);
            else
                photocurrentResponses(iTrial,:,iCone) = pCurrent;
            end
            
%             figure(33);
%             subplot(2,1,1);
%             plot(temporalSupportSeconds, coneExcitation);
%             subplot(2,1,2);
%             plot(temporalSupportSeconds, squeeze(photocurrentResponses(iTrial,:,iCone)));
%             pause
        end
    end
    fprintf('Computed photocurrent response over a time duration of %d msec\n', numel(temporalSupportSeconds)*dt*1000);
%     photocurrentResponses = photocurrentResponses(:,idx,:);
end

function impulseResponse = computePhotocurrentImpulseResponse(timeAxis, cutoffFlag)

    % Temporary solution until we program photocurrent reponse computation
    % in @cMosaic
    load('/Users/qiyuanfeng/Documents/MATLAB/projects/ISETBasicEyeMovements/osLinearFilters25CDM2.mat', 'osLinearFilters', 'dt');
    t = (1:size(osLinearFilters,1))*dt;
    
    if notDefined('cutoffFlag'), cutoffFlag = 'false'; end
    switch cutoffFlag
        case 'true'
            impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), timeAxis);
        case 'false'
            impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), t);
    end
end