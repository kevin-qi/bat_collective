%% Script for the analysis of Voronoi Areas

%close all;

import voronoi_utils.*;
import utils.*;

%Parameters
x2=2.8; x1=-2.8; y2=2.8;  y1=-2.8;  z1=0; z2=2.30;             %Flight volume coordinates
edges_d = {x1:(x2-x1)/10:x2 y1:(y2-y1)/10:y2};                 %Edges for density histogram
bowl = [-0.29, 0.05, 0.45];                                      %x,y,z bowl 1
Fs = 100;                                                      %resampling frequency (Hz) for common time
n_tags = 5;
bat_nms = ['Dai'; 'Den'; 'Dia'; 'Dor'; 'Dum'; 'Ran'];
bat_pairs = nchoosek(1:n_tags,2);   bat_pair_nms = [bat_nms(bat_pairs(:,1),:), '-'.*ones(length(bat_pairs),1), bat_nms(bat_pairs(:,2),:)];
bat_clr = lines(n_tags);
v_th = 0.5;                                                    %Velocity threshold (m/s)
TTL_time_diff = [21; 13; 8; 5; 4];
N = 5;
max_distance = norm([x2-x1 y2-y1]);

bounds = {};
bounds.x1 = x1;
bounds.x2 = x2;
bounds.y1 = y1;
bounds.y2 = y2;

%Custom graded colormap
for i = 1:n_tags
    for j = 1:3
        custom_map(:,j,i) = linspace(1,bat_clr(i,j))';
    end
end

if false
    sessions = ['210222'; '210223'; '210224'; '210225'; '210226'; '210301'; '210302'; '210303']%; '210304'; '210305'];
    session_data = load_session_data(sessions);
    for i=1:length(session_data)
        session_data{i}.x1 = x1;
        session_data{i}.x2 = x2;
        session_data{i}.y1 = y1;
        session_data{i}.y2 = y2;
    end
end
if true
%% Pairwise Distance Distribution
figure;
axes = [];
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    pairwise_dist = pairwise_distance(pos);
    bins = linspace(0,max_distance,20);
    
    index = 1;
    for j = 1:N
        for k = j+1:N
            ax = subplot(length(sessions), N*(N-1)/2, index+(i-1)*N*(N-1)/2);
            axes = [axes ax];
            H = histogram(pairwise_dist{j}(:,k), bins, 'Normalization','probability');
            
            %xline(max_distance, 'LineWidth', 1, 'Color', 'r');
            if i == 1
                title(bat_nms(j,:), bat_nms(k,:));
            end
            if i < length(sessions)
                set(gca,'xtick',[])
            else
                set(gca,'xtick',[1,2,3,4,5,6,7,8]);
            end
            if index > 1
                set(gca,'ytick',[])
            end
            if index == 1
                set(gca, 'ytick',[0.1,0.2,0.3,0.4,0.5]);
                ylabel(sprintf('%s-%s', sessions(i,3:4),sessions(i,5:6)), 'fontweight', 'bold');
            end
            index = index+1;
            
            axis tight;
        end
    end
end
linkaxes(axes, 'xy');
ylim([0,0.5]);
sgtitle('Pairwise Distance Distribution');
end

figure;
axes = [];
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    pos{N+1} = rand_bat_bound_pos(length(pos{1}),bounds);
    pairwise_dist = pairwise_distance(pos);
    bins = linspace(0,max_distance,20);
    
    index = 1;
    for j = 1:N
        for k = j+1:N
            ax = subplot(length(sessions)+1, N*(N-1)/2, index+(i-1)*N*(N-1)/2);
            axes = [axes ax];
            H = histogram(pairwise_dist{j}(:,k), bins, 'Normalization','probability');
            
            %xline(max_distance, 'LineWidth', 1, 'Color', 'r');
            if i == 1
                title(bat_nms(j,:), bat_nms(k,:));
            end
            if i < length(sessions)
                set(gca,'xtick',[])
            else
                set(gca,'xtick',[1,2,3,4,5,6,7,8]);
            end
            if index > 1
                set(gca,'ytick',[])
            end
            if index == 1
                set(gca, 'ytick',[0.1,0.2,0.3,0.4,0.5]);
                ylabel(sprintf('%s-%s', sessions(i,3:4),sessions(i,5:6)), 'fontweight', 'bold');
            end
            index = index+1;
            
            axis tight;
        end
    end
    for j = 1:N*(N-1)/2
        ax = subplot(length(sessions)+1, N*(N-1)/2, j+length(sessions)*N*(N-1)/2);
        axes = [axes ax];
        H = histogram(pairwise_dist{N+1}(:,1), bins, 'Normalization','probability'); 
        if(j == 1)
            ylabel('baseline', 'fontweight', 'bold');
        end
    end
end
linkaxes(axes, 'xy');
ylim([0,0.5]);
sgtitle('Pairwise Distance Distribution With Baseline');

%% Nearest Neighbor Distance Distribution
figure;
axes = [];
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    pairwise_dist = pairwise_distance(pos);
    bins = linspace(0,max_distance,20);
    
    for j = 1:N
        nn_distance = min(pairwise_dist{j}, [], 2);
        ax = subplot(length(sessions), N, j+N*(i-1));
        axes = [axes ax];
        H = histogram(nn_distance, bins, 'Normalization','probability');
        
        xline(max_distance, 'LineWidth', 1, 'Color', 'r');
        if i == 1
            title(bat_nms(j,:));
        end
        if j == 1
            ylabel(sprintf('%s-%s', sessions(i,3:4),sessions(i,5:6)), 'fontweight', 'bold');
        end
    end
    ylim([0,0.5]);
    sgtitle('Nearest Neighbor Distance Distribution');
end
linkaxes(axes, 'xy');


figure;
axes = [];
N=6;
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    pos{N} = rand_bat_bound_pos(length(pos{1}),bounds);
    pairwise_dist = pairwise_distance(pos);
    bins = linspace(0,max_distance,20);
    
    for j = 1:N
        nn_distance = min(pairwise_dist{j}, [], 2);
        ax = subplot(length(sessions), N, j+N*(i-1));
        axes = [axes ax];
        H = histogram(nn_distance, bins, 'Normalization','probability');
        
        xline(max_distance, 'LineWidth', 1, 'Color', 'r');
        if i == 1
            title(bat_nms(j,:));
        end
        if j == 1
            ylabel(sprintf('%s-%s', sessions(i,3:4),sessions(i,5:6)), 'fontweight', 'bold');
        end
    end
    ylim([0,0.5]);
    sgtitle('NN Distance Distribution With Baseline');
end
linkaxes(axes, 'xy');
N=5;


%% Average Distance Distribution
figure;
axes = [];
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    pairwise_dist = pairwise_distance(pos);
    bins = linspace(0,max_distance,20);
    
    for j = 1:N
        avg_distance = mean(pairwise_dist{j}(pairwise_dist{j} ~= 0), 2);
        ax = subplot(length(sessions), N, j+N*(i-1));
        axes = [axes ax];
        H = histogram(avg_distance, bins, 'Normalization','probability');
        
        xline(max_distance, 'LineWidth', 1, 'Color', 'r');
        if i == 1
            title(bat_nms(j,:));
        end
        if j == 1
            ylabel(sprintf('%s-%s', sessions(i,3:4),sessions(i,5:6)), 'fontweight', 'bold');
        end
    end
end
linkaxes(axes, 'xy');
ylim([0,0.5]);
sgtitle('Average Distance Distribution');