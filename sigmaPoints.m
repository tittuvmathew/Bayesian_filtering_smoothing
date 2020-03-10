function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%

    switch type        
        case 'UKF'   
            % your code
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
                
        case 'CKF'
            
            % your code
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
            
        otherwise
            error('Incorrect type of sigma point')
    end

end