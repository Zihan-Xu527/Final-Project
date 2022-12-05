function linear_regression(x, y, yb)

tbl = table(x ,y);
mdl = fitlm(tbl,'linear');
expected_X = 50.6289;
expected_Y = 1.4222;

figure();
plot(mdl)
hold on
scatter(mean(x), mean(y), 'rx')
hold on
scatter(expected_X, mean(yb), 'ro')
hold on
scatter(expected_X, expected_Y, 'ko')
hold off
title('Linear regression')
xlabel('X')
ylabel('Y')