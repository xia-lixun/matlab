close all; clear all; clc;

x = load('tune_positive.dat');
y = load('tune_negative.dat');

z = [x;y];
mu = mean(z);
[coeff,score,latent] = princomp(z);


xhat = bsxfun(@minus, x, mu);
yhat = bsxfun(@minus, y, mu);

xscore = xhat * coeff;
yscore = yhat * coeff;

%%
figure;
plot(xscore(:,4), xscore(:,2), '+'); hold;
plot(yscore(:,4), yscore(:,2), 's');

% the first two are more significant than the second two
%Xapprox = score(:,1:16) * coeff(:,1:16)';
%Xapprox = bsxfun(@plus, mu, Xapprox); % add the mean back in

%%
plot(Xapprox(:,1), Xapprox(:,2), '+')
%% tutorial

Y = randn(100, 2); % 100 observations, 2 features
W = [1 1; 1 -1; 2 1; 2 -1]'; % loading matrix
X = Y * W + 0.1 * randn(100, 4); % add randn noise to simulate experimental error

[coeff, score, latent] = princomp(X);

% score is the matrix of principle components, it should pull out factors
% very close to the original Y variables.
corr(score(:,1),Y(:,1)) % returns -0.9981
corr(score(:,2),Y(:,2)) % returns -0.9830

% score * coeff' will recreate your original data minus its mean,
% i.e. pc x w' = X - mu(X)
mu = mean(X);
xhat = bsxfun(@minus, X, mu);
norm(score * coeff' - xhat)

% coeff is orthogonal, we also have (X-mu)*coeff = score
% so, X = mu + score * coeff'

% To get an approximation to your original data, 
% you can start dropping columns from the computed principal components. 
% To get an idea of which columns to drop, we examine the ev variable
latent

% the first two are more significant than the second two
Xapprox = score(:,1:2) * coeff(:,1:2)';
Xapprox = bsxfun(@plus, mu, Xapprox); % add the mean back in

plot(Xapprox(:,1),X(:,1),'.'); hold on; plot([-4 4],[-4 4])
xlabel('Approximation'); ylabel('Actual value'); grid on;

% a coarser approximation, we just use the first principle component
Xapprox = score(:,1) * coeff(:,1)';
Xapprox = bsxfun(@plus,mu,Xapprox);
plot(Xapprox(:,1),X(:,1),'.'); hold on; plot([-4 4],[-4 4])
xlabel('Approximation'); ylabel('Actual value'); grid on;

% Finally, you might want to see how much of the variance is explained by 
% each of the factors. You can do this using the ev variables:
100*latent/sum(latent)

% So the first component explains 78% of the variance, 
% the next component explains about 22%, and the tiny remainder is 
% explained in the final two components.