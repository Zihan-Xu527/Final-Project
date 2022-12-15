function linear_regression(x, y)

tbl = table(x ,y);
mdl = fitlm(tbl,'linear');

figure()
plot(mdl)
xlabel('X')
ylabel('Y')