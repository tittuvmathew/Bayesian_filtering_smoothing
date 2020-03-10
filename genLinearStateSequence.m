function X = genLinearStateSequence(x_0, P_0, A, Q, N)
%GENLINEARSTATESEQUENCE generates an N-long sequence of states using a 
%    Gaussian prior and a linear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N+1] State vector sequence
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Your code here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialise X 
n = numel(x_0);
X = [];
% draw the initial guess
x0 = mvnrnd(x_0',P_0);          % 1 x n
X = [X, x0'];                   % n x 1 

% Use only the process model 
for k = 2:N+1
   xk = A * X(:,k-1) + (mvnrnd(zeros(1,n),Q))';
   X(:,k) = xk;    
end
end