%% Script for analysing destination choice of flights away from other bats or after foraging (isolating flights)
% This excludes any flights that are towards another bat.

% Imports
import utils.*;
import flight_utils.*;

% Parameters
x2=2.8; x1=-2.8; y2=2.8;  y1=-2.8;  z1=0; z2=2.30;             %Flight volume coordinates
edges_d = {x1:(x2-x1)/10:x2 y1:(y2-y1)/10:y2};                 %Edges for density histogram
bowl = [-0.29, 0.05, 0.45];                                      %x,y,z bowl 1
Fs = 100;                                                      %resampling frequency (Hz) for common time
n_tags = 5;
bat_nms = ['Dai'; 'Den'; 'Dia'; 'Dor'; 'Dum'];
bat_pairs = nchoosek(1:n_tags,2);   bat_pair_nms = [bat_nms(bat_pairs(:,1),:), '-'.*ones(length(bat_pairs),1), bat_nms(bat_pairs(:,2),:)];
bat_clr = lines(n_tags);
v_th = 0.5;                                                    %Velocity threshold (m/s)
TTL_time_diff = [21; 13; 8; 5; 4];
N = 5;

%% Load data
sessions = ['210222'; '210223'; '210224'; '210225'; '210226'; '210301'; '210302'; '210303'; '210304'; '210305'];
sessions = ['210304'; '210305'];
session_data = load_session_data(sessions);
for i=1:length(session_data)
    session_data{i}.x1 = x1;
    session_data{i}.x2 = x2;
    session_data{i}.y1 = y1;
    session_data{i}.y2 = y2;
end

%% Identify isolating flights
for i=1:length(session_data)
    disp(i);
end