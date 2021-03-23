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

verbose = true;

%% Load data
if false
    sessions = ['210222'; '210223'; '210224'; '210225'; '210226'; '210301'; '210302'; '210303'; '210304'; '210305'];
    session_data = load_session_data(sessions);
    for i=1:length(session_data)
        session_data{i}.x1 = x1;
        session_data{i}.x2 = x2;
        session_data{i}.y1 = y1;
        session_data{i}.y2 = y2;
    end
end

%% Identify isolating flights

clear isolating_flights;
figure(1);
hold on;
for i=1:length(session_data)
    flight_data = characterize_flights(session_data{i}, 'bowl', 0.6, 0.6, 0);
    flights{i} = flight_data;
    disp(i);
    for b_num = 1:N
        ind = 1;
        for f_num = 1:length(flight_data{b_num})
            flight = flight_data{b_num}(f_num);
            if((flight.to_bat | flight.to_edge) & (~flight.concurrent & ~flight.self))
                isolating_flights{i, b_num}(ind) = flight;
                ind = ind + 1;
                
                if verbose
                    start_pos = session_data{i}.tag_data_filt{b_num}(flight.start, 3:4);
                    end_pos = session_data{i}.tag_data_filt{b_num}(flight.end, 3:4);
                    
                    scatter(start_pos(1), start_pos(2), 'MarkerFaceColor', 'red');
                    scatter(end_pos(1), end_pos(2), 'MarkerFaceColor', 'blue');
                end
            end
        end
    end
end
disp('done');
hold off;
title('Isolating flights (start at red, end at blue)');
shg;


%% Simulation

% 1. For each flight by a bat (all session and all bats pooled)
% 2. The 'feature vector' is the positions of the non-flying bats, the
% label is the destination of the flying bat. The prediction(s) are the
% destinations (could be multiple valid destinations) of the simulated
% bats. The error is the euclidean distance between predicted and real
% destination.
% 3. Calculate the landscapes for each metric given the configuration of
% the flight.
% 4. Instantiate the different metric optimizing agents
% 5. Get the destinations of a metric-optimal bat


for i=1:length(session_data)
    disp(i);
    for b_num = 1:N
        disp(b_num);
        for f_num = 1:length(isolating_flights{i, b_num})
            flight = isolating_flights{i, b_num}(f_num);
            
            % Indices of the non-flying bats
            stat_b_nums = linspace(1, N, N);
            stat_b_nums = stat_b_nums(stat_b_nums ~= b_num);
            
            % Configuration (feature vector X) of non-flying bats
            X = zeros(2, N);
            for j=1:length(stat_b_nums)
                s_num = stat_b_nums(j);
                X(:, s_num) = session_data{i}.tag_data_filt{s_num}(flight.start, 3:4).';
            end
            
            % Starting position of flying bat
            r0 = session_data{i}.tag_data_filt{b_num}(flight.start, 3:4);
            
            % Ending position of flying bat
            r1 = session_data{i}.tag_data_filt{b_num}(flight.end, 3:4);
            
            trajectory = session_data{i}.tag_data_filt{b_num}(flight.start+1:flight.end-1, 3:4);
            
            
            figure;
            subplot(2,4,1);
            hold on;
            scatter(X(1,stat_b_nums), X(2,stat_b_nums));
            scatter(r0(1), r0(2), 'MarkerFaceColor', 'red');
            scatter(r1(1), r1(2), 'MarkerFaceColor', 'blue');
            scatter(trajectory(:, 1), trajectory(:, 2), 0.03);
            xlim([x1 x2]);
            ylim([x1 x2]);
            axis square;
            hold off;
            shg;
            %assert true == false
            
            bounds.x1 = session_data{i}.x1;
            bounds.x2 = session_data{i}.x2;
            bounds.y1 = session_data{i}.y1;
            bounds.y2 = session_data{i}.y2;
            
            %% Voro
            Z = compute_metric_landscape(X, b_num, 'voro', bounds, 40);
            
            subplot(2,4,2);
            hm = heatmap(linspace(1,40,40),linspace(1,40,40),(Z));
            title('Personal Voronoi Area');
            shg;
           
            top_row = Z(1,:);
            right_col = Z(:,end).';
            bot_row = flip(Z(end,:));
            left_col = flip(Z(:,1)).';
            C = horzcat(top_row, right_col, bot_row, left_col);
            
            [min_dist, argmin] = closest_point_on_perimeter(r1, bounds, 40);
            disp(argmin);
            subplot(2,4,6);
            plot(C);
            hold on;
            x_index = linspace(1,length(C),length(C));
            scatter(x_index(argmin), C(argmin), 'MarkerFaceColor', 'blue');
            [min_dist, argmin] = closest_point_on_perimeter(r0, bounds, 40);
            scatter(x_index(argmin), C(argmin), 'MarkerFaceColor', 'red');
            %% Geo voro
            Z = compute_metric_landscape(X, b_num, 'voro-geo-mean', bounds, 40);
            
            subplot(2,4,3);
            hm = heatmap(linspace(1,40,40),linspace(1,40,40),(Z));
            title('Voronoi Geometric Mean');
            shg;
           
            top_row = Z(1,:);
            right_col = Z(:,end).';
            bot_row = flip(Z(end,:));
            left_col = flip(Z(:,1)).';
            C = horzcat(top_row, right_col, bot_row, left_col);
            
            [min_dist, argmin] = closest_point_on_perimeter(r1, bounds, 40);
            disp(argmin);
            subplot(2,4,7);
            plot(C);
            hold on;
            x_index = linspace(1,length(C),length(C));
            scatter(x_index(argmin), C(argmin), 'MarkerFaceColor', 'blue');
            [min_dist, argmin] = closest_point_on_perimeter(r0, bounds, 40);
            scatter(x_index(argmin), C(argmin), 'MarkerFaceColor', 'red');
            %% NN
            
            Z = compute_metric_landscape(X, b_num, 'nearest-neighbor', bounds, 40);
            
            subplot(2,4,4);
            hm = heatmap(linspace(1,40,40),linspace(1,40,40),(Z));
            title('Nearest Neighbor');
            shg;
           
            top_row = Z(1,:);
            right_col = Z(:,end).';
            bot_row = flip(Z(end,:));
            left_col = flip(Z(:,1)).';
            C = horzcat(top_row, right_col, bot_row, left_col);
            
            [min_dist, argmin] = closest_point_on_perimeter(r1, bounds, 40);
            subplot(2,4,8);
            plot(C);
            hold on;
            x_index = linspace(1,length(C),length(C));
            scatter(x_index(argmin), C(argmin), 'MarkerFaceColor', 'blue');
            [min_dist, argmin] = closest_point_on_perimeter(r0, bounds, 40);
            scatter(x_index(argmin), C(argmin), 'MarkerFaceColor', 'red');
            savefig(sprintf('../public/week_5/flight_sim/%s_%d_%d.fig', sessions(i,:), b_num, f_num));
        end
    end
end
