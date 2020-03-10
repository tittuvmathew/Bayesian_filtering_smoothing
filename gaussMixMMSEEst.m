function [ xHat ] = gaussMixMMSEEst( w, mu, sigma2 )
%GAUSSMIXMMSEEST calculates the MMSE estimate from a Gaussian mixture
%density with multiple components.
%
%Input
%   W           Vector of all the weights
%   MU          Vector containing the means of all components
%   SIGMA2      Vector containing the variances of all components
%
%Output
%   xHat        MMSE estimate

%YOUR CODE HERE
xHat  = 0 ;
nGMM = numel(w);
for i = 1:nGMM
   xHat = xHat + w(i) * mu(i); 
end
end