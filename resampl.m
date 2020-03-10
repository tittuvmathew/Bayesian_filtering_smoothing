function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

% Your code here!
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