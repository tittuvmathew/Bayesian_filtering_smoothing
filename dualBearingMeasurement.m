function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
%DUOBEARINGMEASUREMENT calculates the bearings from two sensors, located in 
%s1 and s2, to the position given by the state vector x. Also returns the
%Jacobian of the model at x.
%
%Input:
%   x           [n x 1] State vector, the two first element are 2D position
%   s1          [2 x 1] Sensor position (2D) for sensor 1
%   s2          [2 x 1] Sensor position (2D) for sensor 2
%
%Output:
%   hx          [2 x 1] measurement vector
%   Hx          [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.
% Your code here
n = size(x,1);
% For Sensor 1 
y1 =  x(2) - s1(2);   % yk - sy
x1 =  x(1) - s1(1); 
hx_1 = atan2(y1,x1);
Hx_1_x = -y1/(x1^2 + y1^2); 
Hx_1_y = x1/(x1^2 + y1^2); 

% For Sensor 2
y2 =  x(2) - s2(2); 
x2 =  x(1) - s2(1); 
hx_2 = atan2(y2,x2); 
Hx_2_x = -y2/(x2^2 + y2^2); 
Hx_2_y = x2/(x2^2 + y2^2); 

hx = [hx_1;hx_2];
Hx = zeros(2,n); 
Hx(1:2,1:2) = [Hx_1_x, Hx_1_y; ...
               Hx_2_x, Hx_2_y];
end