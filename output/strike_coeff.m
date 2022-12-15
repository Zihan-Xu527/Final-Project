
clc; clear all; close all;

% Figure 1.1: strike = 40, n = 50
data = csvread('strike_40.csv');
x = data(:,1);
y = data(:,2);
linear_regression(x,y);

% Figure 1.2: strike = 55, n = 50
data = csvread('strike_55.csv');
x = data(:,1);
y = data(:,2);
linear_regression(x, y);


% Figure 1.3: strike = 70, n = 50
data = csvread('strike_70.csv');
x = data(:,1);
y = data(:,2);
linear_regression(x,y);