function [outputArg1] = kalman_filter(z, a, dt)
%KALMAN_FILTER  Performs a kalman filter on bat positional data.
%   Detailed explanation goes here
%   - z: measured state vector timeseries, state vector := [x_pos, x_vel, y_pos, y_vel]
%   - a: acceleration vector timeseries, accel vector := [x_acc, y_acc]
%   - dt: sampling interval

F = [1 dt 0 0; 
     0 1 0 0; 
     0 0 1 dt; 
     0 0 0 1];
% G = [dt^2/2; dt; dt^2/2; dt];
H = [1 0 0 0; 0 0 1 0]; % Measure only position

% TODO: Calculate sigma_a for x_acc and y_acc separately
sigma_a = mean(sqrt(var(a)));

Q = [dt^4/4 dt^3/2 dt^4/4 dt^3/2;
     dt^3/2 dt^2 dt^3/2 dt^2;
     dt^4/4 dt^3/2 dt^4/4 dt^3/2;
     dt^3/2 dt^2 dt^3/2 dt^2;] * (sigma_a^2);

% TODO: Estimate position measurement variance
sigma_x = 0.2;
sigma_y = 0.2;

R = [sigma_x^2 0; 0 sigma_y^2];

N = length(z);
x0 = [z(1,1) 0 z(1,2) 0].';
P0 = [sigma_x^2 0 0 0;
      0 2*sigma_x 0 0;
      0 0 sigma_y^2 0;
      0 0 0 2*sigma_y];

x_k_k = x0;
P_k_k = P0;
for t=1:N
    % Supposed ground truth evolution (We do not know this ofcourse)
    % x_truth(:,t+1) = F*x_truth(:,t) + G*a(t);
    % z(t) = H*x_truth(:,t) + randn(1,1)*sqrt(R);
    
    
    % Prediction
    x_k1_k = F*x_k_k; % State estimate prediction
    P_k1_k = F*P_k_k*F.' + Q; % State estimate covariance
    
    % Measurement
    z_k = z(t);
    
    % Update
    m_k = z_k - H*x_k1_k; % Innovation
    S_k = H*P_k1_k*H.' + R; % Innovation covariance
    K_k = P_k1_k*H.'*inv(S_k); % Optimal Kalman Gain
    
    x_k1_k1 = x_k1_k + K_k*m_k; % Updated state estimate
    P_k1_k1 = (eye(4,4) - K_k*H)*P_k1_k; % Updated state covariance
    m_k1 = z_k - H*x_k1_k1;
    
    x_estimate(:,t) = x_k1_k1;
    
    x_k_k = x_k1_k1;
    P_k_k = P_k1_k1;
    
end

outputArg1 = x_estimate;
end

