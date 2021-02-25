function [outputArg1] = kalman_filter_bats(x, z, )
%KALMAN_FILTER  Performs a kalman filter on bat positional data.
%   Detailed explanation goes here


for t=1:N
    x_truth(:,t+1) = F*x_truth(:,t) + G*a(t);
    y(t) = H*x_truth(:,t) + randn(1,1)*sqrt(R);
    
    % Prediction
    x_k1_k = F*x_k_k; % State estimate prediction
    P_k1_k = F*P_k_k*F.' + Q; % State estimate covariance
    
    % Measurement
    y_k = y(t);
    
    % Update
    m_k = y_k - H*x_k1_k; % Innovation
    S_k = H*P_k1_k*H.' + R; % Innovation covariance
    K_k = P_k1_k*H.'*inv(S_k); % Optimal Kalman Gain
    
    x_k1_k1 = x_k1_k + K_k*m_k; % Updated state estimate
    P_k1_k1 = (eye(2,2) - K_k*H)*P_k1_k; % Updated state covariance
    m_k1 = y_k - H*x_k1_k1;
    
    x_estimate(:,t) = x_k1_k1;
    
    x_k_k = x_k1_k1;
    P_k_k = P_k1_k1;
    
end

outputArg1 = x_estimate;
end

