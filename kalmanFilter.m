function [X, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%

%% Parameters
N = size(Y,2);      % Number of time steps
n = length(x_0);    % Parametric dimension of state
m = size(Y,1);      % Parametric dimension of measurement

%% Data allocation
X = zeros(n,N);
P = zeros(n,n,N);
% Y, x_0, P_0, A, Q, H, R)
XX = [];
XX = [XX, x_0];  % n x 1 
PP = cell(N+1,1);
PP{1} = P_0;

% Code begins here
for k = 2:N+1
   
   % Prediction 
   x_pred = A * XX(:,k-1);           % Prediction mean 
   P_pred =  A * PP{k-1} * A' + Q;           % Prediction covariance
   
   % Update step
   vk = Y(:,k-1) - H * x_pred;          % Innovation
   Sk = H * P_pred * H' + R;            % Innovation covariance 
   Kk = P_pred * H' * inv(Sk);          % Kalman gain
   
   x_upd = x_pred + Kk * vk;
   P_upd = P_pred - Kk * Sk * Kk';
   
   XX = [XX, x_upd]; 
   PP{k} = P_upd;
   
   X(:,k-1) = x_upd; 
   P(:,:,k-1) = P_upd;    
end



end

