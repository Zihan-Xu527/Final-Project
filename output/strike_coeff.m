
clc; clear all; close all;

data = csvread('strike_40.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
linear_regression(x,y, yb);


data = csvread('strike_50.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
linear_regression(x,y, yb);


data = csvread('strike_60.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
linear_regression(x,y,yb);


data = csvread('strike_70.csv');
x = data(:,1);
y = data(:,2);
yb = data(:,3);
linear_regression(x,y,yb);