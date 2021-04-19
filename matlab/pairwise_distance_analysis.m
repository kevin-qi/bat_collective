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

% Spatial bins for position histograms
bins_x = linspace(x1,x2,20);
bins_y = linspace(y1,y2,20);

%Custom graded colormap
for i = 1:n_tags
    for j = 1:3
        custom_map(:,j,i) = linspace(1,bat_clr(i,j))';
    end
end

%% Flags

% Load Data
load_data_flag = false;

% Calculate pairwise statistics
calc_pairwise_stats_flag = true;


%% Load Data
if load_data_flag
    sessions = ['210222'; '210223'; '210224'; '210225'; '210226'; '210301'; '210302'; '210303'; '210304'; '210305'; '210308'; '210309'; '210310';'210311';'210315';'210316';'210317';'210318';'210319'];
    session_data = load_session_data(sessions);
    for i=1:length(session_data)
        session_data{i}.x1 = x1;
        session_data{i}.x2 = x2;
        session_data{i}.y1 = y1;
        session_data{i}.y2 = y2;
    end
end

%% Calculate pairwise statistics
if calc_pairwise_stats_flag

N_bins = 20;
bins = linspace(0,max_distance,N_bins);
base_counts = zeros(length(sessions), N*(N-1)/2, length(bins)-1);
true_counts = zeros(length(sessions), N*(N-1)/2, length(bins)-1);
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    for j = 1:N
        rand_pos{j} = pos{j}(randperm(length(pos{j})),:);
    end
    pos{N+1} = rand_bat_bound_pos(length(pos{1}),bounds);
    pairwise_dist{i} = pairwise_distance(pos);
    shuffled_pairwise_dist{i} = pairwise_distance(rand_pos);

    index = 1;
    for j = 1:N
        for k = j+1:N

            base_counts(i,index,:) = histcounts(shuffled_pairwise_dist{i}{j}(:,k), bins);
            true_counts(i,index,:) = histcounts(pairwise_dist{i}{j}(:,k), bins);

            base_counts(i,index,:) = base_counts(i,index,:) / sum(base_counts(i,index,:));
            true_counts(i,index,:) = true_counts(i,index,:) / sum(true_counts(i,index,:));

            diff_counts(i,index,:) = true_counts(i,index,:) - base_counts(i,index,:);

            [base_ecdf{i, index}.f base_ecdf{i, index}.x] = ecdf(shuffled_pairwise_dist{i}{j}(:,k));
            [true_ecdf{i, index}.f true_ecdf{i, index}.x] = ecdf(pairwise_dist{i}{j}(:,k));
            
            % Get CDF(x=1)
            [min_val, argmin] = min(abs(true_ecdf{i, index}.x - 1));
            true_proximity(i, index) = true_ecdf{i,index}.f(argmin);
            
            [min_val, argmin] = min(abs(base_ecdf{i, index}.x - 1));
            base_proximity(i, index) = base_ecdf{i,index}.f(argmin);
            
            proximity_metric(i, index) = true_proximity(i, index) - base_proximity(i, index);

            counts_j = histcounts2(pos{j}(:,2), pos{j}(:,1), bins_y, bins_x);
            counts_j = counts_j / sum(sum(counts_j));

            counts_k = histcounts2(pos{k}(:,2), pos{k}(:,1), bins_y, bins_x);
            counts_k = counts_k / sum(sum(counts_k));

            bin_centers_x = bins_x + diff(bins_x(1:2))/2;
            bin_centers_x = bin_centers_x(1:end-1);
            bin_centers_y = bins_y + diff(bins_y(1:2))/2;
            bin_centers_y = bin_centers_y(1:end-1);
            grid_coords_x = meshgrid(bin_centers_x, bin_centers_y);
            grid_coords_y = flip(meshgrid(bin_centers_x, bin_centers_y).',1);
            j_peaks{i,index} = [grid_coords_x(counts_j>0.1) grid_coords_y(counts_j>0.1)].';
            k_peaks{i,index} = [grid_coords_x(counts_k>0.1) grid_coords_y(counts_k>0.1)].';
            index = index + 1;
        end
    end
end

save for_angelo/pairwise_distances.mat pairwise_dist shuffled_pairwise_dist
save for_angelo/pairwise_dist_counts.mat base_counts true_counts diff_counts
save for_angelo/empirical_CDF.mat base_ecdf true_ecdf
save for_angelo/proximity_measures.mat true_proximity base_proximity proximity_metric
end