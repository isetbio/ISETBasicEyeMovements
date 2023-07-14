function visualizeClassificationMatrices(classificationMatrix, classificationMatrixProjection, taskIntervals)
    figure(); clf;
    subplot(2,1,1);
    visualizeClassificationMatrix(classificationMatrix, taskIntervals, 'cone index');
    subplot(2,1,2);
    visualizeClassificationMatrix(classificationMatrixProjection, [], 'principal component no');
end
