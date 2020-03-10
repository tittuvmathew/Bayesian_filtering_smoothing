function [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type)
%NONLINEARKALMANFILTER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence for times 1,...,N
%   x_0         [n x 1] Prior mean for time 0
%   P_0         [n x n] Prior covariance
%   f                   Motion model function handle
%                       [fx,Fx]=f(x) 
%                       Takes as input x (state) 
%                       Returns fx and Fx, motion model and Jacobian evaluated at x
%   Q           [n x n] Process noise covariance
%   h                   Measurement model function handle
%                       [hx,Hx]=h(x,T) 
%                       Takes as input x (state), 
%                       Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%

% Your code here. If you have good code for the Kalman filter, you should re-use it here as
% much as possible.
m = size(Y,1); 
N = size(Y,2); 
n = size(x_0,1); 
% initialising 
xf = zeros(n,N);    Pf = zeros(n,n,N);       % Posterior mean and covariance 
xp = zeros(n,N);    Pp = zeros(n,n,N);       % Predicted mean and covariance

switch type
    case 'EKF'       
        % Loop over the time step 
        for i = 1:N
           if i == 1    % First step 
              % EKF prediction
              [F,G] = f(x_0); 
              xpred = F; 
              Ppred = G * P_0 * G' + Q;  
              % Store the predicted mean and covariance for EKF
              xp(:,i) = xpred; 
              Pp(:,:,i) = Ppred; 
              
              % EKF update 
              [hx,Hx] = h(xpred); 
              Sk = Hx * Ppred * Hx' + R; 
              Kk = Ppred * Hx' * inv(Sk); 
              xpost = xpred + Kk * (Y(:,i) - hx); 
              Ppost = Ppred - Kk * Sk * Kk';   
              % Store the updated mean and covariance for EKF (posterior)
              xf(:,i) = xpost;
              Pf(:,:,i) = Ppost;  
           else
              % EKF prediction 
              [F,G] = f(xf(:,i-1));
              xpred = F; 
              Ppred = G * Pf(:,:,i-1) * G' + Q;  
              % Store the predicted mean and covariance for EKF
              xp(:,i) = xpred; 
              Pp(:,:,i) = Ppred; 
              
              % EKF update 
              [hx,Hx] = h(xpred); 
              Sk = Hx * Ppred * Hx' + R; 
              Kk = Ppred * Hx' * inv(Sk); 
              xpost = xpred + Kk * (Y(:,i) - hx); 
              Ppost = Ppred - Kk * Sk * Kk';   
              % Store the updated mean and covariance for EKF (posterior)
              xf(:,i) = xpost;
              Pf(:,:,i) = Ppost;                 
           end
        end
    % Case regarding UKF         
    case 'UKF'
        % Loop over the time step 
        for i = 1 : N 
            if i == 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%   UKF prediction
                n = size(x_0,1); 
                W0 = 1 - (n/3); 
                W = [W0 , ((1-W0)/(2*n))*ones(1,2*n)]; 
                SP = zeros(n,2*n+1); 
                SP(:,1) = x_0; 
                Psqrt = sqrtm(P_0);       % matrix square root of P 
                % Loop over each column of P for plus signs 
                for j = 1:n
                    SP(:,j+1) = x_0 + sqrt(n/(1-W0)) * Psqrt(:,j);
                end
                % Loop over each column of P for minus signs 
                for j = 1:n
                    SP(:,n+1+j) = x_0 - sqrt(n/(1-W0)) * Psqrt(:,j); 
                end
                % Compute the approximated predicted mean and covariance using UKF
                xpred = zeros(n,1); 
                for j = 1 : 2*n+1
                    xpred = xpred + W(j) * f(SP(:,j)) ;
                end
                Ppred = zeros(n,n);
                for j = 1:2*n+1
                    Ppred = Ppred + W(j) * (f(SP(:,j)) - xpred) * (f(SP(:,j)) - xpred)' ;
                end
                Ppred = Ppred + Q; 
                xp(:,i) = xpred; 
                Pp(:,:,i) = Ppred;
                
                %%%%%%%%%%%%%%%%%%%%%%%   UKF update / posterior 
                n = size(x_0,1); 
                W0 = 1 - (n/3); 
                W = [W0 , ((1-W0)/(2*n))*ones(1,2*n)]; 
                SP = zeros(n,2*n+1); 
                SP(:,1) = xpred; 
                Psqrt = sqrtm(Ppred);       % matrix square root of P 
                % Loop over each column of P for plus signs 
                for j = 1:n
                    SP(:,j+1) = xpred + sqrt(n/(1-W0)) * Psqrt(:,j);
                end
                % Loop over each column of P for minus signs 
                for j = 1:n
                    SP(:,n+1+j) = xpred - sqrt(n/(1-W0)) * Psqrt(:,j); 
                end
                m = size(Y,1);
                yhat = zeros(m,1); 
                for j = 1 : 2*n+1
                    yhat = yhat + W(j) * h(SP(:,j)); 
                end
                Pxy = zeros(n,m); 
                for j = 1:2*n+1
                    Pxy = Pxy + W(j) * ((SP(:,j) - xpred) * (h(SP(:,j)) - yhat)') ; 
                end
                Sk = zeros(m,m);
                for j = 1:2*n+1
                    Sk = Sk + W(j) *  (h(SP(:,j)) - yhat) * (h(SP(:,j)) - yhat)';
                end
                Sk = Sk + R; 
                % Approximated mean and covariance 
                xpost = xpred + Pxy * inv(Sk) * (Y(:,i) - yhat); 
                Ppost = Ppred - Pxy * inv(Sk) * Pxy';   
                xf(:,i) = xpost;
                Pf(:,:,i) = Ppost;  
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%   UKF prediction
                n = size(x_0,1);
                W0 = 1 - (n/3); 
                W = [W0 , ((1-W0)/(2*n))*ones(1,2*n)]; 
                SP = zeros(n,2*n+1); 
                
                SP(:,1) = xf(:,i-1); 
                Psqrt = sqrtm(Pf(:,:,i-1));        % matrix square root of P 
                
                % Loop over each column of P for plus signs 
                for j = 1:n
                    SP(:,j+1) = xf(:,i-1) + sqrt(n/(1-W0)) * Psqrt(:,j);
                end
                % Loop over each column of P for minus signs 
                for j = 1:n
                    SP(:,n+1+j) = xf(:,i-1) - sqrt(n/(1-W0)) * Psqrt(:,j); 
                end
                % Compute the approximated predicted mean and covariance using UKF
                xpred = zeros(n,1); 
                for j = 1 : 2*n+1
                    xpred = xpred + W(j) * f(SP(:,j)) ;
                end
                Ppred = zeros(n,n);
                for j = 1:2*n+1
                    Ppred = Ppred + W(j) * (f(SP(:,j)) - xpred) * (f(SP(:,j)) - xpred)' ;
                end
                Ppred = Ppred + Q; 
                xp(:,i) = xpred; 
                Pp(:,:,i) = Ppred;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%% UKF update / posterior 
                W0 = 1 - (n/3); 
                W = [W0 , ((1-W0)/(2*n))*ones(1,2*n)]; 
                SP = zeros(n,2*n+1); 
                SP(:,1) = xpred; 
                Psqrt = sqrtm(Ppred);       % matrix square root of P 
                % Loop over each column of P for plus signs 
                for j = 1:n
                    SP(:,j+1) = xpred + sqrt(n/(1-W0)) * Psqrt(:,j);
                end
                % Loop over each column of P for minus signs 
                for j = 1:n
                    SP(:,n+1+j) = xpred - sqrt(n/(1-W0)) * Psqrt(:,j); 
                end
                
                yhat = zeros(m,1); 
                for j = 1 : 2*n+1
                    yhat = yhat + W(j) * h(SP(:,j)); 
                end
                Pxy = zeros(n,m); 
                for j = 1:2*n+1
                    Pxy = Pxy + W(j) * ((SP(:,j) - xpred) * (h(SP(:,j)) - yhat)') ; 
                end
                Sk = zeros(m,m);
                for j = 1:2*n+1
                    Sk = Sk + W(j) *  (h(SP(:,j)) - yhat) * (h(SP(:,j)) - yhat)';
                end
                Sk = Sk + R; 
                % Approximated mean and covariance of the posterior 
                xpost = xpred + Pxy * inv(Sk) * (Y(:,i) - yhat); 
                Ppost = Ppred - Pxy * inv(Sk) * Pxy';   
                xf(:,i) = xpost;
                Pf(:,:,i) = Ppost; 
            end
        end   
    case 'CKF' 
        % Loop over the time step 
        for i = 1 : N 
            if i == 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%   CKF prediction
                W = [(1/(2*n))*ones(1,2*n)];
                Psqrt = sqrtm(P_0);       % matrix square root of P 
                SP = zeros(n,2*n); 
                % Loop over each column of P for plus signs 
                for j = 1:n
                    SP(:,j) = x_0 + sqrt(n) * Psqrt(:,j);
                end
                % Loop over each column of P for minus signs 
                for j = 1:n
                    SP(:,n+j) = x_0 - sqrt(n) * Psqrt(:,j); 
                end
                % Compute the approximated predicted mean and covariance
                xpred = zeros(n,1); 
                for j = 1 : 2*n
                    xpred = xpred + W(j) * f(SP(:,j)) ;
                end
                Ppred = zeros(n,n);
                for j = 1:2*n
                    Ppred = Ppred + W(j) * (f(SP(:,j)) - xpred) * (f(SP(:,j)) - xpred)' ;
                end
                Ppred = Ppred + Q; 
                xp(:,i) = xpred; 
                Pp(:,:,i) = Ppred;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CKF Update/ Posterior 
                W = [(1/(2*n))*ones(1,2*n)];
                Psqrt = sqrtm(Ppred);       % matrix square root of P 
                SP = zeros(n,2*n); 
                % Loop over each column of P for plus signs 
                for j = 1:n
                    SP(:,j) = xpred + sqrt(n) * Psqrt(:,j);
                end
                % Loop over each column of P for minus signs 
                for j = 1:n
                    SP(:,n+j) = xpred - sqrt(n) * Psqrt(:,j); 
                end
                % Update mean and covariance 
                yhat = zeros(m,1); 
                for j = 1 : 2*n
                    yhat = yhat + W(j) * h(SP(:,j)); 
                end 
                Pxy = zeros(n,m); 
                for j = 1:2*n
                    Pxy = Pxy + W(j) * ((SP(:,j) - xpred) * (h(SP(:,j)) - yhat)') ; 
                end
                Sk = zeros(m,m);
                for j = 1:2*n
                    Sk = Sk + W(j) *  (h(SP(:,j)) - yhat) * (h(SP(:,j)) - yhat)';
                end
                Sk = Sk + R; 
                % Approximated mean and covariance 
                xpost = xpred + Pxy * inv(Sk) * (Y(:,i) - yhat); 
                Ppost = Ppred - Pxy * inv(Sk) * Pxy';   
                xf(:,i) = xpost;
                Pf(:,:,i) = Ppost; 
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%   CKF prediction
                W = [(1/(2*n))*ones(1,2*n)];
                Psqrt = sqrtm(Pf(:,:,i-1));         % matrix square root of P 
                SP = zeros(n,2*n); 
                % Loop over each column of P for plus signs 
                for j = 1:n
                    SP(:,j) = xf(:,i-1) + sqrt(n) * Psqrt(:,j);
                end
                % Loop over each column of P for minus signs 
                for j = 1:n
                    SP(:,n+j) = xf(:,i-1) - sqrt(n) * Psqrt(:,j); 
                end
                % Compute the approximated predicted mean and covariance
                xpred = zeros(n,1); 
                for j = 1 : 2*n
                    xpred = xpred + W(j) * f(SP(:,j)) ;
                end
                Ppred = zeros(n,n);
                for j = 1:2*n
                    Ppred = Ppred + W(j) * (f(SP(:,j)) - xpred) * (f(SP(:,j)) - xpred)' ;
                end
                Ppred = Ppred + Q; 
                xp(:,i) = xpred; 
                Pp(:,:,i) = Ppred;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CKF Update/ Posterior 
                W = [(1/(2*n))*ones(1,2*n)];
                Psqrt = sqrtm(Ppred);       % matrix square root of P 
                SP = zeros(n,2*n); 
                % Loop over each column of P for plus signs 
                for j = 1:n
                    SP(:,j) = xpred + sqrt(n) * Psqrt(:,j);
                end
                % Loop over each column of P for minus signs 
                for j = 1:n
                    SP(:,n+j) = xpred - sqrt(n) * Psqrt(:,j); 
                end
                % Update mean and covariance 
                yhat = zeros(m,1); 
                for j = 1 : 2*n
                    yhat = yhat + W(j) * h(SP(:,j)); 
                end 
                Pxy = zeros(n,m); 
                for j = 1:2*n
                    Pxy = Pxy + W(j) * ((SP(:,j) - xpred) * (h(SP(:,j)) - yhat)') ; 
                end
                Sk = zeros(m,m);
                for j = 1:2*n
                    Sk = Sk + W(j) *  (h(SP(:,j)) - yhat) * (h(SP(:,j)) - yhat)';
                end
                Sk = Sk + R; 
                % Approximated mean and covariance 
                xpost = xpred + Pxy * inv(Sk) * (Y(:,i) - yhat); 
                Ppost = Ppred - Pxy * inv(Sk) * Pxy';   
                xf(:,i) = xpost;
                Pf(:,:,i) = Ppost; 
            end
        end       
    otherwise
         error('Incorrect type of non-linear Kalman filter')
end
end