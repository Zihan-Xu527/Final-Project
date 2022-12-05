clc; clear all; close all

expected_X = 50.6289;
expected_Y = 1.4222;

data1 = csvread('sample_size_10.csv');
x1 = data1(:,1);
y1 = data1(:,2);
yb1 = data1(:,3);
tbl1 = table(x1 , y1);
mdl1 = fitlm(tbl1,'linear');

figure;
plot(mdl1)
hold on
scatter(mean(x1), mean(y1), 'rx')
hold on
scatter(expected_X, mean(yb1), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off
title('n = 10')
xlabel('X')
ylabel('Y')


data2 = csvread('sample_size_50.csv');
x2 = data2(:,1);
y2 = data2(:,2);
yb2 = data2(:,3);
tbl2 = table(x2 , y2);
mdl2 = fitlm(tbl2,'linear');

figure();
plot(mdl2)
hold on
scatter(mean(x2), mean(y2), 'rx')
hold on
% scatter(x2, yb2, 'o')
% hold on
scatter(expected_X, mean(yb2), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off
title('n = 50')
xlabel('X')
ylabel('Y')

data3 = csvread('sample_size_100.csv');
x3 = data3(:,1);
y3 = data3(:,2);
yb3 = data3(:,3);
tbl3 = table(x3 , y3);
mdl3 = fitlm(tbl3,'linear');

figure();
plot(mdl3)
hold on
scatter(mean(x3), mean(y3), 'rx')
hold on
% scatter(x3, yb3, 'o')
% hold on
scatter(expected_X, mean(yb3), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off
title('n = 100')
xlabel('X')
ylabel('Y')

data4 = csvread('sample_size_1000.csv');
x4 = data4(:,1);
y4 = data4(:,2);
yb4 = data4(:,3);
tbl4 = table(x4 , y4);
mdl4 = fitlm(tbl4,'linear');

figure();
plot(mdl4)
hold on
scatter(mean(x4), mean(y4), 'rx')
hold on
% scatter(x4, yb4, 'o')
% hold on
scatter(expected_X, mean(yb4), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off
title('n = 1000')
xlabel('X')
ylabel('Y')