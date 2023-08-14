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