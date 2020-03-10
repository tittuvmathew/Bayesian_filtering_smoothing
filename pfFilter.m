function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

% Your code here, please. 
% If you want to be a bit fancy, then only store and output the particles if the function
% is called with more than 2 output arguments.
n = size(x_0,1);
m = size(Y,1); 
K = size(Y,2); 

Xp = zeros(n,N,K);
Wp = zeros(N,K);
xfp = zeros(n,K);
Pfp = zeros(n,n,K);

% Generate N samples from prior 
X_0_N = mvnrnd(x_0',P_0,N);   %  N x n 
W_0_N = (1./N) * ones(1,N);   % All with equal weights 

for i = 1:K
   % if its the first step 
   %i
   if i == 1
      if bResample == 0     % no resample 
          [X_k, W_k] = pfFilterStep(X_0_N', W_0_N, Y(:,i) , proc_f, proc_Q, meas_h, meas_R); 
          Xp(:,:,i) = X_k; 
          Wp(:,i) = W_k';    
          % Compute posterior mean and coavriance of particles 
          xfp(:,i) = Xp(:,:,i) * Wp(:,i);
          Pfp(:,:,i) = (Xp(:,:,i) - xfp(:,i)) * ((Xp(:,:,i) - xfp(:,i))'.* Wp(:,i));  
      elseif bResample == 1     % resample 
          %i
          % First resample and then perform prediction 
          [X_k, W_k, j] = resampl_new(X_0_N', W_0_N);  
          [X_k, W_k] = pfFilterStep(X_k, W_k, Y(:,i) , proc_f, proc_Q, meas_h, meas_R); 
          
          Xp(:,:,i) = X_k;
          Wp(:,i) = W_k';
          % Compute posterior mean and coavriance of particles
          xfp(:,i) = Xp(:,:,i) * Wp(:,i);
          Pfp(:,:,i) = (Xp(:,:,i) - xfp(:,i)) * ((Xp(:,:,i) - xfp(:,i))'.* Wp(:,i));
      end     
   elseif i > 1
      if bResample == 0 
          Wpi = (Wp(:,i-1))';
          [X_k, W_k] = pfFilterStep(Xp(:,:,i-1),Wpi, Y(:,i),proc_f,proc_Q,meas_h,meas_R); 
          Xp(:,:,i) = X_k; 
          Wp(:,i) = W_k';
          % Compute posterior mean of particle filter 
          xfp(:,i) = Xp(:,:,i) * Wp(:,i);
          Pfp(:,:,i) = (Xp(:,:,i) - xfp(:,i)) * ((Xp(:,:,i) - xfp(:,i))'.* Wp(:,i)); 
      elseif bResample == 1 
          %i
          Wpi = (Wp(:,i-1))';
          % First resample and then perform prediction 
          [X_k, W_k, j] = resampl_new(Xp(:,:,i-1),Wpi);  
          [X_k, W_k] = pfFilterStep(X_k,W_k,Y(:,i),proc_f,proc_Q,meas_h,meas_R); 
          Xp(:,:,i) = X_k; 
          Wp(:,i) = W_k';
          xfp(:,i) = Xp(:,:,i) * Wp(:,i);
          Pfp(:,:,i) = (Xp(:,:,i) - xfp(:,i)) * ((Xp(:,:,i) - xfp(:,i))'.* Wp(:,i));         
      end
   end  
end
end

function [Xr, Wr, j] = resampl_new(X, W)
    n = size(X,1);      % dimension of a sample 
    N = size(X,2);      % number of samples 
    Wdiv = [];
    Wdiv = [Wdiv 0];
    Wi = 0; 
    % Loop over the weights 
    indx_remov = [];
    for i = 1:N
        Wi = Wi + W(i); 
        if W(i) == 0 
            indx_remov = [indx_remov , i];
        end      
        Wdiv = [Wdiv Wi];
    end
    % Take the unique values 
    Wdiv = unique(Wdiv); 
    indx_kept = setdiff(1:N,indx_remov);   % Samples having non-zero initial weights arranged in sorted order
    % uniform number generation 
    Rnumb = rand(1,N); 
    %Rnumb = [0.65,0.03,0.84,0.93];
    counter = zeros(1,numel(Wdiv)-1);
    % Loop over each uniform samples and check in which bin they fall 
    for i = 1:numel(Rnumb)
        samplei = Rnumb(i);
        % Loop over the sub-division 
        indxi = [];
        for j = 1:numel(Wdiv)-1
            if Wdiv(j) < samplei &  samplei < Wdiv(j+1)  % Check whether we need to add <= and >=
                indxi = [indxi, j]; 
            end
        end
        % increment the counter 
        atemp = zeros(1,numel(Wdiv)-1);
        atemp(indxi) = 1; 
        counter = counter + atemp; 
    end
    % Equal weights after resampling 
    W_per_sample = 1/N;
    Wr = ones(1,N) * W_per_sample;
    % Generate Xr 
    Xr = [];
    j = [];
    % Loop over the counter 
    for i = 1:numel(Wdiv)-1
        if counter(i) > 0 
            Xr = [Xr , repmat(X(:,indx_kept(i)),1,counter(i))];
            j = [j , repmat(indx_kept(i),1,counter(i))];
        end
    end
end

function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
    % Copy your code from previous task!
    X_k_T = mvnrnd((proc_f(X_kmin1))',proc_Q);   %  N x n 
    X_k = X_k_T';   %   n x N 
    W_k = W_kmin1.*(mvnpdf(yk',(meas_h(X_k))',meas_R))'; 
    W_k = W_k./sum(W_k);
end