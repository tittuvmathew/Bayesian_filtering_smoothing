function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
%PFFILTERSTEP Compute one filter step of a SIS/SIR particle filter.
%
% Input:
%   X_kmin1     [n x N] Particles for state x in time k-1
%   W_kmin1     [1 x N] Weights for state x in time k-1
%   y_k         [m x 1] Measurement vector for time k
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%
% Output:
%   X_k         [n x N] Particles for state x in time k
%   W_k         [1 x N] Weights for state x in time k

% Your code here!

n = size(X_kmin1,1); 
N = size(X_kmin1,2); 

% Loop over the samples 
X_k = zeros(n,N);  W_k = zeros(1,N);

% Loop over the samples
%for i = 1:N
%    % Get the i'th sample 
%    xi = X_kmin1(:,i); 
%    xk_i = mvnrnd((proc_f(xi))',proc_Q);  % New sample from sampling density 
%    X_k(:,i) = xk_i';    
%    wnewi = W_kmin1(i) * mvnpdf(yk',(meas_h(xk_i'))',meas_R); 
%    W_k(1,i) = wnewi;
%end
%wsum = sum(W_k); 
%W_k = W_k ./ wsum;
    X_k_T = mvnrnd((proc_f(X_kmin1))',proc_Q);   %  N x n 
    X_k = X_k_T';   %   n x N 
    W_k = W_kmin1'.*mvnpdf(yk',(meas_h(X_k))',meas_R); 
    W_k = W_k'./sum(W_k);
end