load('/Users/qiyuanfeng/Documents/MATLAB/projects/ISETBasicEyeMovements/osLinearFilters25CDM2.mat');
t = (1:size(osLinearFilters,1))*0.00033; % 0.0033
figure()
plot(t, osLinearFilters(:,1),'b-','LineWidth',2)
hold on;
plot(t, osLinearFilters(:,2),'LineWidth',1.9, 'color', "#77AC30")
plot(t, osLinearFilters(:,3),'r-','LineWidth',1.8)
yticks(-0.02:0.04:0.18);
xlabel('Time (seconds)');
xlim([0 0.12])
ax = gca; % current axes
ax.FontSize = 12;
ax.LineWidth = 1.5;
% ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% ax.XLim = [0 0.12];
hold off
box off
