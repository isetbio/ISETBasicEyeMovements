%% IRF
irf = impulseResponse;

figure();
axis = 1:400;
plot(axis, irf);

%% Conv
x = zeros(1, 400);
x(100:200) = 1;

y = conv(x, irf, 'full');
y = y(1:length(x));

figure();
plot(axis, x); hold on;
plot(axis, y); 

