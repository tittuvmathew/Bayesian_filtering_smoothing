function Y = genLinearMeasurementSequence(X, H, R)
%GENLINEARMEASUREMENTSEQUENCE generates a sequence of observations of the state 
% sequence X using a linear measurement model. Measurement noise is assumed to be 
% zero mean and Gaussian.
%
%Input:
%   X           [n x N+1] State vector sequence. The k:th state vector is X(:,k+1)
%   H           [m x n] Measurement matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% your code here
N = size(X,2) - 1; 
n = size(X,1);
m = size(H,1);
Y = [];
for k = 1:N
    Yk = H * X(:,k+1) + (mvnrnd(zeros(1,m),R))'; 
    Y = [Y, Yk];
end   
    
end