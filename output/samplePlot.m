clc; clear all; close all;

expected_X = 50.6289;
expected_Y = 1.4222;

data = csvread('sample_size_10.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
linear_regression(x, y, yb)
hold on
scatter(mean(x), mean(y), 'rx')
hold on
scatter(expected_X, mean(yb), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off


data = csvread('sample_size_50.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
linear_regression(x, y, yb)
hold on
scatter(mean(x), mean(y), 'rx')
hold on
scatter(expected_X, mean(yb), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off

data = csvread('sample_size_100.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
linear_regression(x, y, yb)
hold on
scatter(mean(x), mean(y), 'rx')
hold on
scatter(expected_X, mean(yb), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off

data = csvread('sample_size_1000.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
linear_regression(x, y, yb)
hold on
scatter(mean(x), mean(y), 'rx')
hold on
scatter(expected_X, mean(yb), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off





