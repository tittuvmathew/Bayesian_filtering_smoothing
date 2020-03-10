function X = genNonLinearStateSequence(x_0, P_0, f, Q, N)
%GENLINEARSTATESEQUENCE generates an N+1-long sequence of states using a 
%    Gaussian prior and a linear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N+1] State vector sequence
%

% Your code here

% Sample from the prior first 
x_1 = mvnrnd(x_0',P_0); 
n = size(x_0,1); 
X = zeros(n,N+1);
X(:,1) = x_1';

for k = 2:N+1
   xtemp = f(X(:,k-1)) + (mvnrnd(zeros(1,n),Q))';
   X(:,k) = xtemp;
end


end