function [curvefit, gof] = fitPsychometric(x,y)
    [curvefit, gof] = fit(x(:), y(:), 'a+b/(1+exp(-c*(x-d)))');
    iteration = 1;
    while gof.rsquare < 0.9 && iteration < 5
        fprintf('poor curve fitting, trying again. Iteration: %d\n', iteration);
        [curvefit, gof] = fit(x(:), y(:), 'a+b/(1+exp(-c*(x-d)))');
        iteration = iteration + 1;
    end
    if gof.rsquare < 0.9
        fprintf('Warning: poor auto curve fitting after 5 iterations\n');
    end
    
    sse = gof.sse;
    rsquare = gof.rsquare;
    setThreshold = 75;
    find = @(x) setThreshold - curvefit(x);
    threshold = fzero(find, 0);
    
    fprintf('Detection Threshold is %.4f with a r-square of %.4f\n', threshold, rsquare);
    
    figure()
    plot(x,y,'.', 'MarkerSize', 12);
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