%%Parameters
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
N_bins = 20;
bins = linspace(0,max_distance,N_bins);
sessions = ['210222'; '210223'; '210224'; '210225'; '210226'; '210301'; '210302'; '210303'; '210304'; '210305'; '210308'; '210309'; '210310';'210311';'210315';'210316';'210317';'210318';'210319'];

%% Load data
if false
% Contains pairwise distances and shuffled pairwise distances
load('../data/processed/pairwise_analysis/pairwise_distances.mat');

% Contains histogram bin counts (normalized) for true, shuffled, and
% difference histograms
load('../data/processed/pairwise_analysis/pairwise_dist_counts.mat');

% Contains the empirical cumulative distribution of true and shuffled
load('../data/processed/pairwise_analysis/empirical_CDF.mat');

% Contains proximity measures (fraction of data with dist < 1) of true and
% shuffled data as well as their difference (the proximity metric)
load('../data/processed/pairwise_analysis/proximity_measures.mat');

% Load session data
load('../data/processed/session_data.mat');

end

%% Interbat distance
figure;
subset_ind = [1 3 5 7 9 11 13 15 17 19];
t = tiledlayout(length(subset_ind),N*(N-1)/2, ...
            'TileSpacing','Compact', ...
            'Padding','Compact');

for i = 1:length(sessions) % i-th session
    index = 1; % Index ranges from 1 to 10, denotes the bat pair.
    for j = 1:N % j-th bat
        for k = j+1:N % k-th bat (not equal to j-th bat)
            if ismember(i, subset_ind)
                nexttile

                % Plotting the true, observed distribution
                bar(reshape(true_counts(i,index,:), [1 length(bins)-1]), 'BarWidth', 1, 'FaceColor', 'k');
                ylim([0,0.6]);

                % label x-ticks every 2 meters
                xticks([0, N_bins*2/max_distance, N_bins*4/max_distance, N_bins*6/max_distance, N_bins*8/max_distance]);
                xticklabels([0,2,4,6,8]);

                index = index+1; 
                if i == 1
                    title(sprintf('%s-%s', bat_nms(j,:), bat_nms(k,:)));
                end
            end
        end
    end
end
title(t, 'Observed Interbat Distance Distributions');
xlabel(t, 'Distance (m)');
ylabel(t, 'Session');
savefig('../data/processed/pairwise_analysis/results/observed_interbat_distance_distributions_skip.fig');
%saveas(gcf, 'for_angelo/observed_interbat_distance_distributions.png');

%% Interbat distance with overlay
figure;
bins_x = linspace(x1,x2,10);
bins_y = linspace(y1,y2,10);
subset_ind = [1 3 5 7 9 11 13 15 17 19];
t = tiledlayout(length(subset_ind),N*(N-1)/2, ...
            'TileSpacing','Compact', ...
            'Padding','Compact');

for i = 1:length(sessions) % i-th session
    index = 1; % Index ranges from 1 to 10, denotes the bat pair.
    data = session_data{i};
    pos = extract_position(data);
    for j = 1:N % j-th bat
        for k = j+1:N % k-th bat (not equal to j-th bat)
            if ismember(i, subset_ind)
                nexttile

                % Plotting the true, observed distribution
                bar(reshape(true_counts(i,index,:), [1 length(bins)-1]), 'BarWidth', 1, 'FaceColor', 'k');
                ylim([0,0.6]);
                
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
                        xline(norm(j_p - k_p)*19/max_distance, 'LineWidth', 1, 'Color', 'r');
                    end
                end
                
                % label x-ticks every 2 meters
                xticks([0, N_bins*2/max_distance, N_bins*4/max_distance, N_bins*6/max_distance, N_bins*8/max_distance]);
                xticklabels([0,2,4,6,8]);

                index = index+1; 
                if i == 1
                    title(sprintf('%s-%s', bat_nms(j,:), bat_nms(k,:)));
                end
            end
        end
    end
end
title(t, 'Observed Interbat Distance Distributions');
xlabel(t, 'Distance (m)');
ylabel(t, 'Session');
savefig('../data/processed/pairwise_analysis/results/observed_interbat_distance_distributions_overlay_skip.fig');
%saveas(gcf, 'for_angelo/observed_interbat_distance_distributions.png');

%% Interbat shuffled distance
figure;
t = tiledlayout(length(subset_ind),N*(N-1)/2, ...
            'TileSpacing','Compact', ...
            'Padding','Compact');

for i = 1:length(sessions)
    index = 1; % Index ranges from 1 to 10, denotes the bat pair.
    for j = 1:N % j-th bat
        for k = j+1:N % k-th bat (not equal to j-th bat)
            if ismember(i, subset_ind)
                nexttile

                % Plotting the shuffled baseline distribution
                bar(squeeze(mean(base_counts(i,index,:,:),4)), 'BarWidth', 1, 'FaceColor', 'k');
                ylim([0,0.6]);

                % label x-ticks every 2 meters
                xticks([0, N_bins*2/max_distance, N_bins*4/max_distance, N_bins*6/max_distance, N_bins*8/max_distance]);
                xticklabels([0,2,4,6,8]);
               
                index = index+1;
                if i == 1
                    title(sprintf('%s-%s', bat_nms(j,:), bat_nms(k,:)));
                end
            end
        end
    end
end
title(t, 'Baseline Interbat Distance Distributions');
xlabel(t, 'Distance (m)');
ylabel(t, 'Session');
savefig('../data/processed/pairwise_analysis/results/shuffled_interbat_distance_distributions_skip.fig');
%saveas(gcf, 'for_angelo/shuffled_interbat_distance_distributions.png');

%% Interbat distance difference distribution
figure;
t = tiledlayout(length(subset_ind),N*(N-1)/2, ...
            'TileSpacing','Compact', ...
            'Padding','Compact');

for i = 1:length(sessions)
    index = 1;
    for j = 1:N
        for k = j+1:N
            if ismember(i, subset_ind)
                nexttile

                % Plotting the difference (true - baseline) distribution
                bar(squeeze(diff_counts(i,index,:)), 'BarWidth', 1, 'FaceColor', 'k');
                [H, p, ci] = ttest(squeeze(boot_counts(i,index,1,:)),mean(squeeze(base_counts(i,index,1,:))), 'Alpha', 0.01);   
                hold on;
                if(H == 0)
                    alpha(0.1);
                end
                %errorbar( 1+mean(diff(bins))/2, squeeze(diff_counts(i,index,1)), diff(ci)/2);
                ylim([-0.25,0.25]);

                xticks([0, N_bins*2/max_distance, N_bins*4/max_distance, N_bins*6/max_distance, N_bins*8/max_distance]);
                xticklabels([0,2,4,6,8]);
                hold off
                index = index+1;
                if i == 1
                    title(sprintf('%s-%s', bat_nms(j,:), bat_nms(k,:)));
                end
            end
        end
    end
end
title(t, 'Interbat Distance Difference Distributions');
xlabel(t, 'Distance (m)');
ylabel(t, 'Session');
savefig('../data/processed/pairwise_analysis/results/difference_distributions_skip.fig');
%saveas(gcf, 'for_angelo/difference_distributions.png');

if false % skip
%% K-S test true and shuffled pairwise distances
% K-S test between pairwise_dist and shuffled_pairwise_dist
for i = 1:length(sessions)
    index = 1;
    for j = 1:N
        for k = j+1:N
            [h(i,index),p(i,index)] = kstest2(pairwise_dist{i}{j}(:,k), shuffled_pairwise_dist{i}{j}(:,k), 'Alpha', 0.001);
            index = index + 1;
        end
    end
end
end



%% Observed & shuffled interbat distance ECDF
H(1) = figure;
t = tiledlayout(length(sessions),N*(N-1)/2, ...
            'TileSpacing','Compact', ...
            'Padding','Compact');

for i = 1:length(sessions)
    index = 1;
    for j = 1:N
        for k = j+1:N
            nexttile
            
            % Plot of the observed ECDF
            stairs(true_ecdf{i,index}.x, true_ecdf{i,index}.f);
            hold on
            % Plot of the baseline ECDF
            stairs(base_ecdf{i,index}.x, base_ecdf{i,index}.f);
            hold off;
            
            % K-S test between observed and shuffled
            %[h,p] = kstest2(pairwise_dist{i}{j}(:,k), shuffled_pairwise_dist{i}{j}(:,k));
            %title(sprintf('%d',h));
            %ylim([0,1]);
            
            %xticks([0, N_bins*2/max_distance, N_bins*4/max_distance, N_bins*6/max_distance, N_bins*8/max_distance]);
            %xticklabels([0,2,4,6,8]);
            % Critical value separating 'proximal' and 'distal' fractions
            % of the distribution
            xline(1, 'LineWidth', 1, 'Color', 'r');
            
            index = index+1;
            if i == 1
                title(sprintf('%s-%s', bat_nms(j,:), bat_nms(k,:)));
            end
        end
    end
end
title(t, 'Interbat Distance Cumulative Distributions');
xlabel(t, 'Distance (m)');
ylabel(t, 'Session');
lgd = legend('Observed', 'Baseline');
lgd.Layout.Tile = 'East';
%saveas(gcf, 'for_angelo/empirical_CDF.png');
savefig(H, '../data/processed/pairwise_analysis/results/emprical_CDF.fig', 'compact')

%% Proximity metric
H(1) = figure;
t = tiledlayout(N*(N-1)/2,1, ...
            'TileSpacing','Compact', ...
            'Padding','Compact');

index = 1;
for j = 1:N
    for k = j+1:N
        nexttile

        % Plot of the observed ECDF
        plot(proximity_metric(:,index));
        ylim([-0.2 0.2]);
        xlim([1 19]);
        % K-S test between observed and shuffled
        %[h,p] = kstest2(pairwise_dist{i}{j}(:,k), shuffled_pairwise_dist{i}{j}(:,k));
        %title(sprintf('%d',h));
        %ylim([0,1]);
        yline(0);
        %xticks([0, N_bins*2/max_distance, N_bins*4/max_distance, N_bins*6/max_distance, N_bins*8/max_distance]);
        %xticklabels([0,2,4,6,8]);
        % Critical value separating 'proximal' and 'distal' fractions
        % of the distribution
        %xline(1, 'LineWidth', 1, 'Color', 'r');

        index = index+1;
        if i == 1
            title(sprintf('%s-%s', bat_nms(j,:), bat_nms(k,:)));
        end
    end
end
sgtitle(t, 'Proximity metric (Observed proximity - shuffled proximity)');
xlabel(t, 'Distance (m)');
ylabel(t, 'Session');
%saveas(gcf, 'for_angelo/proximity_metric.png');
savefig('../data/processed/pairwise_analysis/results/proximity_metric.fig');


%% Network Graphs
H(1) = figure;
t = tiledlayout(4,5, ...
            'TileSpacing','Compact', ...
            'Padding','Compact');
for i = 1:length(sessions)
    nexttile
    index = 1;
    A = zeros(N);
    for j = 1:N
        for k = j+1:N
            A(j,k) = proximity_metric(i,index);
            index = index + 1;
        end
    end
    A = A + A.'; % Graph adjacency matrix
    G = graph(A);   
    G.Edges.Weight = squeeze(proximity_metric(i,:)).';
    G.Edges.LWidths = 50*abs(G.Edges.Weight); %7*G.Edges.Weight/max(G.Edges.Weight);\
    edge_colors = zeros(N*(N-1)/2,3);
    edge_colors(G.Edges.Weight > 0, :) = repmat([0 0.4470 0.7410], size(edge_colors(G.Edges.Weight > 0, :),1),1);
    edge_colors(G.Edges.Weight <= 0, :) = repmat([0.8500 0.3250 0.0980], size(edge_colors(G.Edges.Weight <= 0, :),1),1);
    
    p_NG = plot(G, 'EdgeColor', edge_colors, 'NodeLabel', num2cell(bat_nms(1:N,:), 2));
    p_NG.LineWidth = G.Edges.LWidths;  p_NG.MarkerSize = 7;  %p_NG.NodeLabelColor = bat_clr;
    p_NG.NodeFontSize = 15;    p_NG.NodeFontWeight = 'bold';
    title(sessions(i,:));
    
end
sgtitle(t, 'Proximity metric Network Graphs');
%saveas(gcf, 'for_angelo/proximity_metric.png');
savefig('../data/processed/pairwise_analysis/results/proximity_metric_network_graphs.fig');
figure();

