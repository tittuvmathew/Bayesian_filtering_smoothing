function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Predicted mean
%   P           [n x n] Predicted covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%

    switch type
        case 'EKF'            
            % Your EKF update here
            [hx,Hx]=h(x); 
            Sk = Hx * P * Hx' + R; 
            Kk = P * Hx' * inv(Sk); 
            x = x + Kk * (y - hx); 
            P = P - Kk * Sk * Kk';              
        case 'UKF'    
            % Your UKF update here
            n = size(x,1); 
            W0 = 1 - (n/3); 
            W = [W0 , ((1-W0)/(2*n))*ones(1,2*n)]; 
            SP = zeros(n,2*n+1); 
            SP(:,1) = x; 
            Psqrt = sqrtm(P);       % matrix square root of P 
            % Loop over each column of P for plus signs 
            for i = 1:n
                SP(:,i+1) = x + sqrt(n/(1-W0)) * Psqrt(:,i);
            end
            % Loop over each column of P for minus signs 
            for i = 1:n
               SP(:,n+1+i) = x - sqrt(n/(1-W0)) * Psqrt(:,i); 
            end
            % Update mean and covariance 
            m = size(y,1); 
            yhat = zeros(m,1); 
            for i = 1 : 2*n+1
               yhat = yhat + W(i) * h(SP(:,i)); 
            end
            Pxy = zeros(n,m); 
            for i = 1:2*n+1
               Pxy = Pxy + W(i) * ((SP(:,i) - x) * (h(SP(:,i)) - yhat)') ; 
            end
            Sk = zeros(m,m);
            for i = 1:2*n+1
               Sk = Sk + W(i) *  (h(SP(:,i)) - yhat) * (h(SP(:,i)) - yhat)';
            end
            Sk = Sk + R; 
            % Approximated mean and covariance 
            x = x + Pxy * inv(Sk) * (y - yhat); 
            P = P - Pxy * inv(Sk) * Pxy';           
        case 'CKF'    
            % Your CKF update here
            n = size(x,1); 
            W = [(1/(2*n))*ones(1,2*n)];
            Psqrt = sqrtm(P);       % matrix square root of P 
            SP = zeros(n,2*n); 
            % Loop over each column of P for plus signs 
            for i = 1:n
                SP(:,i) = x + sqrt(n) * Psqrt(:,i);
            end
            % Loop over each column of P for minus signs 
            for i = 1:n
               SP(:,n+i) = x - sqrt(n) * Psqrt(:,i); 
            end
            % Update mean and covariance 
            m = size(y,1); 
            yhat = zeros(m,1); 
            for i = 1 : 2*n
               yhat = yhat + W(i) * h(SP(:,i)); 
            end
            Pxy = zeros(n,m); 
            for i = 1:2*n
               Pxy = Pxy + W(i) * ((SP(:,i) - x) * (h(SP(:,i)) - yhat)') ; 
            end
            Sk = zeros(m,m);
            for i = 1:2*n
               Sk = Sk + W(i) *  (h(SP(:,i)) - yhat) * (h(SP(:,i)) - yhat)';
            end
            Sk = Sk + R; 
            % Approximated mean and covariance 
            x = x + Pxy * inv(Sk) * (y - yhat); 
            P = P - Pxy * inv(Sk) * Pxy';            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end

