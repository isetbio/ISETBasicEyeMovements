%% Target psychometric curve
axis = 0 : 0.01 : 100;
p = wblcdf(axis, 50, 5);

figure(1);
plot(axis, p, 'k', 'LineWidth', 2);
xlabel('Stimulus');
ylabel('Prob');

%% Simulate response 
stim = rand(1, 100) * 100;
prob = wblcdf(stim, 50, 5);

choice = zeros(size(prob));
for idx = 1:length(prob)
    z = rand();
    choice(idx) = z <= prob(idx);
end

hold on;
scatter(stim, choice);

%% MLE estimate
objective = @(para) - wbl_llhd(para(1), para(2), stim, choice);
[esti, fval] = fmincon(objective, [50, 5], [], [], [], [], [1, 1], []);

%% Plot
p_est = wblcdf(axis, esti(1), esti(2));
plot(axis, p_est, 'r--', 'LineWidth', 2);