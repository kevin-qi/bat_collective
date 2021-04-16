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
    sessions = ['210222'; '210223'; '210224'; '210225'; '210226'; '210301'; '210302'; '210303'; '210304'; '210305'; '210308'; '210309'; '210310';'210311';'210315';'210316';'210317';'210318';'210319'];
    session_data = load_session_data(sessions);
    for i=1:length(session_data)
        session_data{i}.x1 = x1;
        session_data{i}.x2 = x2;
        session_data{i}.y1 = y1;
        session_data{i}.y2 = y2;
    end
end

if true
%% Individual Bat Positions
figure;
axes = [];
for i = 1:round(length(sessions)/2 -1)
    data = session_data{2*i};
    pos = extract_position(data);
    bins_x = linspace(x1,x2,10);
    bins_y = linspace(y1,y2,10);
    
    index = 1;
    for j = 1:N
        ax = subplot(round(length(sessions)/2 -1), N, (i-1)*N + j);
        axes = [axes ax];
        counts = histcounts2(pos{j}(:,2), pos{j}(:,1), bins_y, bins_x);
        counts = counts / sum(sum(counts));
        
        imagesc(bins_y, bins_x, counts);
        axis equal;
        xlim([x1 x2]);
        ylim([y1 y2]);
        caxis([0 1]);
        colormap hot;
        colorbar();
        
        if i == 1
            title(bat_nms(j,:));
        end
        if index == 1
            ylabel(sprintf('%s-%s', sessions(2*i,3:4),sessions(2*i,5:6)), 'fontweight', 'bold');
        end
        index = index+1;
    end
end
linkaxes(axes, 'xy');
sgtitle('Individual Bat Position Distribution');


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
                set(gca,'xtick',[1,2,3,4,5,6,7,8]);
            else
                set(gca,'xtick',[1,2,3,4,5,6,7,8]);
            end
            if index > 1
                set(gca, 'ytick',[0.1,0.2,0.3,0.4,0.5]);
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
    for j = 1:N
        rand_pos{j} = pos{j}(randperm(length(pos{j})),:);
    end
    pos{N+1} = rand_bat_bound_pos(length(pos{1}),bounds);
    pairwise_dist = pairwise_distance(pos);
    rand_pairwise_dist = pairwise_distance(rand_pos);
    bins = linspace(0,max_distance,20);
    
    index = 1;
    for j = 1:N
        for k = j+1:N
            ax = subplot(length(sessions)+1, N*(N-1)/2, index+(i-1)*N*(N-1)/2);
            axes = [axes ax];
            H = histogram(rand_pairwise_dist{j}(:,k), bins, 'Normalization','probability', 'FaceColor', 'k');

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
            j_peaks = [grid_coords_x(counts_j>0.1) grid_coords_y(counts_j>0.1)].';
            k_peaks = [grid_coords_x(counts_k>0.1) grid_coords_y(counts_k>0.1)].';
            
            peak_distances = [];
            for j_p = j_peaks
                for k_p = k_peaks
                    peak_distances = [peak_distances norm(j_p - k_p)];
                    xline(norm(j_p - k_p), 'LineWidth', 1, 'Color', 'r');
                end
            end
            
            
                    
            %xline(max_distance, 'LineWidth', 1, 'Color', 'r');
            if i == 1
                title(bat_nms(j,:), bat_nms(k,:));
            end
            if i < length(sessions)
                set(gca,'xtick',[1,2,3,4,5,6,7,8]);
            else
                set(gca,'xtick',[1,2,3,4,5,6,7,8]);
            end
            if index > 1
                set(gca, 'ytick',[0.1,0.2,0.3,0.4,0.5]);
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
sgtitle('Randomized Pairwise Distance Distribution With Baseline');

figure;
axes = [];
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    for j = 1:N
        rand_pos{j} = pos{j}(randperm(length(pos{j})),:);
    end
    pos{N+1} = rand_bat_bound_pos(length(pos{1}),bounds);
    pairwise_dist = pairwise_distance(pos);
    rand_pairwise_dist = pairwise_distance(rand_pos);
    bins = linspace(0,max_distance,20);
    
    index = 1;
    for j = 1:N
        for k = j+1:N
            ax = subplot(length(sessions)+1, N*(N-1)/2, index+(i-1)*N*(N-1)/2);
            axes = [axes ax];
            H = histogram(pairwise_dist{j}(:,k), bins, 'Normalization','probability', 'FaceColor', 'k');
            %hold on;
            %H2 = histogram(rand_pairwise_dist{j}(:,k), bins, 'Normalization','probability', 'FaceColor', 'k', 'FaceAlpha', 1);
            %hold off;
            %xline(max_distance, 'LineWidth', 1, 'Color', 'r');
            if i == 1
                title(bat_nms(j,:), bat_nms(k,:));
            end
            if i < length(sessions)
                set(gca,'xtick',[1,2,3,4,5,6,7,8]);
            else
                set(gca,'xtick',[1,2,3,4,5,6,7,8]);
            end
            if index > 1
                set(gca, 'ytick',[0.1,0.2,0.3,0.4,0.5]);
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
    % (5, N_samples)
    for j = 1:N
        ax = subplot(length(sessions), N, j+N*(i-1));
        axes = [axes ax];
        [nn_distance, argmin] = min(pairwise_dist{j}, [], 2);
        session_counts = zeros(length(bins)-1, N);
        for k = 1:N
            [counts, edges] = histcounts(nn_distance(argmin == k), bins);
            session_counts(:, k) = counts.';
        end
        
        bar(session_counts/sum(sum(session_counts, 2)), 'stacked');
        %H = histogram(nn_distance, bins, 'Normalization','probability');
        
        xticks([0, 20*2/max_distance, 20*4/max_distance, 20*6/max_distance, 20*8/max_distance]);
        xticklabels([0,2,4,6,8]);
        if i == 1
            title(bat_nms(j,:));
        end
        if i == 1 & j == 5
            legend('Dai', 'Den', 'Dia', 'Dor', 'Dum');
        end
        if j == 1
            ylabel(sprintf('%s-%s', sessions(i,3:4),sessions(i,5:6)), 'fontweight', 'bold');
        end
    end
    ylim([0,0.5]);
    sgtitle('Nearest Neighbor Distance Distribution');
end
linkaxes(axes, 'xy');



%% Nearest Neighbor First Bin
figure;
axes = [];
session_counts = zeros(N, N, length(sessions));
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    pairwise_dist = pairwise_distance(pos);
    bins = linspace(0,max_distance,10);
    % (5, N_samples)
    for j = 1:N
        [nn_distance, argmin] = min(pairwise_dist{j}, [], 2);
        for k = 1:N
            [counts, edges] = histcounts(nn_distance(argmin == k), bins);
            session_counts(j, k, i) = counts(1);
        end
    end
end
for j = 1:N
    ax = subplot(N, 1, j);
    axes = [axes ax];
    area(reshape(session_counts(j,:,:), [N length(sessions)]).');
    ylabel(bat_nms(j,:));
    
    xticks(linspace(1, length(sessions), length(sessions)));
    if j == 1
        legend('Dai', 'Den', 'Dia', 'Dor', 'Dum');
    end
    if j == N
        xlabel('Session');
    end
end
sgtitle('NN Distance Cumulative Distribution < 0.5 meters');

%% Nearest Neighbor First Bin Normalized
figure;
axes = [];
session_counts = zeros(N, N, length(sessions));
for i = 1:length(sessions)
    data = session_data{i};
    pos = extract_position(data);
    pairwise_dist = pairwise_distance(pos);
    bins = linspace(0,max_distance,10);
    % (5, N_samples)
    for j = 1:N
        [nn_distance, argmin] = min(pairwise_dist{j}, [], 2);
        for k = 1:N
            [counts, edges] = histcounts(nn_distance(argmin == k), bins);
            session_counts(j, k, i) = counts(1);
        end
        session_counts(j,:,i) = session_counts(j,:,i)/sum(session_counts(j,:,i));
    end
end
disp(session_counts)
for j = 1:N
    ax = subplot(N, 1, j);
    axes = [axes ax];
    area(reshape(session_counts(j,:,:), [N length(sessions)]).');
    ylabel(bat_nms(j,:));
    xticks(linspace(1, length(sessions), length(sessions)));
    ylim([0 1]);
    if j == 1
        legend('Dai', 'Den', 'Dia', 'Dor', 'Dum');
    end
    if j == N
        xlabel('Session');
    end
end
sgtitle('NN Distance Cumulative Distribution < 0.5 meters (Normalized)');


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


axes = [];
for i = 1:1
    figure;
    data = session_data{i};
    pos = extract_position(data);
    pairwise_dist = pairwise_distance(pos);
    bins = linspace(0,max_distance,20);
    
    index = 1;
    for j = 1:N
        for k = 1:N
            ax = subplot(N, N, (j-1)*N+k);
            axes = [axes ax];
            
            if(j ~= k)
                [counts_j, y_edges, x_edges, binY_j, binX_j] = histcounts2(pos{j}(:,2), pos{j}(:,1), bins_y, bins_x);
                [counts_k, y_edges, x_edges, binY_k, binX_k] = histcounts2(pos{k}(:,2), pos{k}(:,1), bins_y, bins_x);

                mean_dist_j = zeros(length(bins_x)-1);
                mean_dist_k = zeros(length(bins_x)-1);
                for b_x = 1:length(bins_x)-1
                    for b_y = 1:length(bins_y)-1
                        j_in_bin_j = pos{j}(binX_j == b_x & binY_j == b_y, :); % Position of bat j when bat j is in bin b_x
                        k_in_bin_j = pos{k}(binX_j == b_x & binY_j == b_y, :); % Position of bat k when bat j is in bin b_x

                        j_in_bin_k = pos{j}(binX_k == b_x & binY_k == b_y, :); % Position of bat j when bat k is in bin b_x
                        k_in_bin_k = pos{k}(binX_k == b_x & binY_k == b_y, :); % Position of bat k when bat k is in bin b_x

                        mean_dist_j(b_y, b_x) = mean(sqrt(sum((j_in_bin_j - k_in_bin_j).^2,2)));
                        mean_dist_k(b_y, b_x) = mean(sqrt(sum((j_in_bin_k - k_in_bin_k).^2,2)));
                    end
                end
                
                counts_j = histcounts2(pos{j}(:,2), pos{j}(:,1), bins_y, bins_x);
                counts_j = counts_j / sum(sum(counts_j));
                
                counts_k = histcounts2(pos{k}(:,2), pos{k}(:,1), bins_y, bins_x);
                counts_k = counts_k / sum(sum(counts_k));
                
                mask = (counts_j + counts_k) > 0.05;
                inv_mask = ones(size(mask)) - mask;
                inv_mask = inv_mask * 9999;
                
                
                imagesc(bins_y, bins_x, mean_dist_k.*-1);
                axis equal;
                xlim([x1 x2]);
                ylim([y1 y2]);
                caxis([-5 0]);
                colormap parula;
                colorbar();
            else
                counts = histcounts2(pos{j}(:,2), pos{j}(:,1), bins_y, bins_x);
                counts = counts / sum(sum(counts));

                imagesc(bins_y, bins_x, counts);
                axis equal;
                xlim([x1 x2]);
                ylim([y1 y2]);
                caxis([0 1]);
                colormap parula;
                colorbar();
            end
            
        end
    end
end
linkaxes(axes, 'xy');
sgtitle('Conditional Pairwise Distance');

if false
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
end