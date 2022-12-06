clc; clear all; close all;
expected_X = 50.6289;
% expected_Y = 10.6624; %strike = 40
expected_Y = 1.4222;  % strike = 55
% expected_Y = 0.0481762; % strike = 70
data = csvread('sample_size_10000.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
figure()
[counts1, binCenters1] = hist(y, 500);
[counts2, binCenters2] = hist(yb, 500);
plot(binCenters1, counts1, 'g-');
hold on;
plot(binCenters2, counts2, 'b-');
hold on;
line([expected_Y, expected_Y], ylim, 'LineWidth', 2, 'Color', 'r');
grid on;
legend('w/o control variates', 'w/ control variates', 'exact solution')



[f1,y] = ksdensity(y);
[f2,yb] = ksdensity(yb);
figure()
plot(y, f1, 'g-');
hold on
plot(yb,f2, 'b-');
hold on
line([expected_Y, expected_Y], ylim, 'LineWidth', 2, 'Color', 'r');
legend('w/o control variates', 'w/ control variates', 'exact solution')
title('Comparison of distributions')
