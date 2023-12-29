function y = poisson(lambda, numsamples)
sigma = sqrt(1/lambda/2);
x = (sigma*randn(numsamples,1)).^2 + (sigma*randn(numsamples,1)).^2;
y = cumsum(x,1);
end