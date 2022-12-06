function linear_regression(x, y, yb)

tbl = table(x ,y);
mdl = fitlm(tbl,'linear');
expected_X = 50.6289;
expected_Y = 1.4222;

figure()
plot(mdl)
xlabel('X')
ylabel('Y')