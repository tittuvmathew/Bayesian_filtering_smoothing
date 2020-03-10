function [mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N)
%approxGaussianTransform takes a Gaussian density and a transformation 
%function and calculates the mean and covariance of the transformed density.
%
%Inputs
%   MU_X        [m x 1] Expected value of x.
%   SIGMA_X     [m x m] Covariance of x.
%   F           [Function handle] Function which maps a [m x 1] dimensional
%               vector into another vector of size [n x 1].
%   N           Number of samples to draw. Default = 5000.
%
%Output
%   MU_Y        [n x 1] Approximated mean of y.
%   SIGMA_Y     [n x n] Approximated covariance of y.
%   ys          [n x N] Samples propagated through f


if nargin < 4
    N = 5000;
end

%Your code here
xsamples = mvnrnd(mu_x,Sigma_x,N);  % xamples of size N x m   

y_s = [];
for i =1:N
   y_s = [y_s , f(xsamples(i,:)')]; 
end

% Approximated mean of Y
mu_y = mean(y_s,2); 

% Approximated covariance of Y
Sigma_y = zeros(size(mu_y,1),size(mu_y,1)); 
for i = 1: N
   Sigma_y = Sigma_y + (y_s(:,i) - mu_y) * (y_s(:,i) - mu_y)';
end

Sigma_y = Sigma_y ./ (N-1); 

%cov_y = cell(N,1); 
%for i = 1:N
%   cov_y{i} =  (y_s(:,i)-mu_y) * (y_s(:,i)-mu_y)';
%end

%for i = 1:size(mu_y,1)
%   for j = 1:size(mu_y,1)
%      sum = 0; 
%      for k = 1:N
%         sum = sum + cov_y{k}(i,j); 
%      end
%      Sigma_y(i,j) = sum/N; 
%   end
%end
    
end