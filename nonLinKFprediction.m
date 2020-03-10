function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

    switch type
        case 'EKF'
            
            % Your EKF code here
            [F,G] = f(x); 
            x = F; 
            P = G * P * G' + Q; 
            
        case 'UKF'
            % Your UKF code here
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
            % Compute the approximated predicted mean and covariance
            xnew = zeros(n,1); 
            for i = 1 : 2*n+1
               xnew = xnew + W(i) * f(SP(:,i)) ;
            end
            x = xnew; 
            Pnew = zeros(n,n);
            for i = 1:2*n+1
                Pnew = Pnew + W(i) * (f(SP(:,i)) - x) * (f(SP(:,i)) - x)' ;
            end
            Pnew = Pnew + Q; 
            P = Pnew ; 
            
        case 'CKF'            
            % Your CKF code here
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
            % Compute the approximated predicted mean and covariance
            xnew = zeros(n,1); 
            for i = 1 : 2*n
               xnew = xnew + W(i) * f(SP(:,i)) ;
            end
            x = xnew; 
            Pnew = zeros(n,n);
            for i = 1:2*n
                Pnew = Pnew + W(i) * (f(SP(:,i)) - x) * (f(SP(:,i)) - x)' ;
            end
            Pnew = Pnew + Q; 
            P = Pnew ;       
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end