%% Script for the analysis of Voronoi Areas

import voronoi_utils.*;
import utils.*;

%Parameters
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

%Custom graded colormap
for i = 1:n_tags
    for j = 1:3
        custom_map(:,j,i) = linspace(1,bat_clr(i,j))';
    end
end
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

if false
v_areas = {};
stat_durations = {};
for sess_index = 1:length(sessions)
    sess_name = sessions(sess_index, :);
    data = session_data{sess_index};

    data.x1=x1;
    data.x2=x2;
    data.y1=y1;
    data.y2=y2;

    [v_areas{sess_index}, stat_durations{sess_index}] = compute_voronoi_session(data, true, false, false);
end
end

figure;
for i = 1:length(v_areas)
    bat_areas = v_areas{i};
    bins = linspace(0,30,15);
    for j = 1:N
        subplot(length(v_areas),N,j+N*(i-1));
        
        H = histogram(bat_areas(j,:), bins, 'Normalization','pdf');
        xline(31.36/5, 'LineWidth', 1, 'Color', 'r');
        ylim([0 0.3]);
        if i == 1
            title(bat_nms(j,:));
        end
        if j == 1
            ylabel(sprintf('%s-%s', sessions(i,3:4),sessions(i,5:6)), 'fontweight', 'bold');
        end
    end
    
    sgtitle('Voronoi Cell Area Distributions (unique)');
    
end

figure;
for i = 1:length(v_areas) % session
    bat_areas = v_areas{i};
    bins = linspace(0,30,15);
    for j = 1:N % bat
        subplot(length(v_areas),N,j+N*(i-1));
        disp(size(bat_areas(j,:)));
        disp(size(stat_durations{i}));
        
        weighted_areas = []
        for k = 1:length(bat_areas(j,:))
            weighted_areas = [weighted_areas ones(1,stat_durations{i}(k))*bat_areas(j,k)];
        end
        disp(size(weighted_areas));
        disp(sum(stat_durations{i}));
        
        H = histogram(weighted_areas, bins, 'Normalization','pdf');
        xline(31.36/5, 'LineWidth', 1, 'Color', 'r');
        ylim([0 0.3]);
        if i == 1
            title(bat_nms(j,:));
        end
        if j == 1
            ylabel(sprintf('%s-%s', sessions(i,3:4),sessions(i,5:6)), 'fontweight', 'bold');
        end
        
    end
    
    sgtitle('Voronoi Cell Area Distributions (weighted)');
    
end


%% Compute theoretical maximal possible area per bat configurations
if(false)
    N_iter = 10000;
    N_points = 10;
    X = linspace(data.x1, data.x2, N_points);
    Y = linspace(data.y1, data.y2, N_points);

    [A,B] = meshgrid(X, [data.y1 data.y2]);
    [C,D] = meshgrid([data.x1 data.x2], Y);
    R = cat(1,[A(:) B(:)],[C(:) D(:)]);
    R = unique(R, 'rows');
    figure;
    scatter(R(:,1), R(:,2)); % Monte Carlo Position Sampling
    shg;
    area=zeros(length(V), N_iter);

    bat_configuration = {};
    for i=1:N_iter
        random_index = randperm(length(R), 5);

        P = R(random_index, :).';
        bat_configuration{i} = P;
        [A, b, V] = voronoiPolyhedronsCCW(P, [-2.8 -2.8], [2.8 2.8]);
        for j=1:length(V)
            % X & Y coordinates of the j-th voronoi cell vertices
            V_x = V{j}(1,:);
            V_y = V{j}(2,:);

            area(j,i) = polyarea(V_x, V_y);
        end
    end

    figure;
    bins = linspace(0,30,30);
    for i=1:N
        subplot(1,N,i);
        histogram(area(i,:), bins, 'Normalization','probability');
    end

    figure;
    [max_prod, max_ind] = max(prod(area,1));
    disp(max_ind);
    scatter(bat_configuration{max_ind}(1,:), bat_configuration{max_ind}(2,:));
    hold on;
    [A, b, V] = voronoiPolyhedronsCCW(bat_configuration{max_ind}, [-2.8 -2.8], [2.8 2.8]);
    for j=1:length(V)
        % X & Y coordinates of the j-th voronoi cell vertices
        V_x = V{j}(1,:);
        V_y = V{j}(2,:);

        a = polyarea(V_x, V_y);

        % Compute centroid
        cent_x = mean(V_x);
        cent_y = mean(V_y);

        % Complete the polygon by appending the first entry (for
        % plotting purposes)
        V_x(end+1) = V_x(1); 
        V_y(end+1) = V_y(1);

        plot(V_x, V_y);
        text(cent_x, cent_y+0.1, string(a));
    end
    title(string(sum(area(:, max_ind))));
    hold off;
    shg;
end


%% Voronoi Area Contour Plots
N_points = 40;
X = linspace(data.x1, data.x2, N_points);
Y = linspace(data.y1, data.y2, N_points);
[A,B] = meshgrid(X, Y);
XY = [A(:) B(:)];

for i = 1:N
    bat_x_pos(:, i) = data.tag_data_filt{i}(1:350000,3);
    bat_y_pos(:, i) = data.tag_data_filt{i}(1:350000,4);
end


if false
figure;
subplot(3,3,5);

t=180000;
for i = 1:N
    scatter(bat_x_pos(t,i), bat_y_pos(t,i));
    
    dx = -0.05; dy = 0.2; % displacement so the text does not overlay the data points
    text(bat_x_pos(t,i)+dx, bat_y_pos(t,i)+dy, string(i));
    hold on
end
ylim([data.y1, data.y2]);
xlim([data.x1, data.x2]);
hold off

% Personal Voronoi Area Heatmap (bat 1)
subplot(3,3,1);
Z = zeros(N_points);
for i=1:N_points
    for k=1:N_points
        P = cat(1, [X(i) bat_x_pos(t, 2:end)], [Y(k) bat_y_pos(t, 2:end)]);
        [A, b, V] = voronoiPolyhedronsCCW(P, [-2.8 -2.8], [2.8 2.8]);
        
        V_x = V{1}(1,:);
        V_y = V{1}(2,:);

        Z(i,k) = polyarea(V_x, V_y);
    end
end

hm = heatmap(X,Y,Z);
title('Individual Voronoi Area of Bat 1');
hm.YDisplayData=flip(hm.YDisplayData);

XLabels = linspace(data.x1, data.x2, N_points);
CustomXLabels = string(XLabels);
CustomXLabels(mod(XLabels, 1) ~= 0) = " ";
CustomXLabels(1) = string(data.x1);
CustomXLabels(end) = string(data.x2);

YLabels = linspace(data.y1, data.y2, N_points);
CustomYLabels = string(YLabels);
CustomYLabels(mod(YLabels, 1) ~= 0) = " ";
CustomYLabels(1) = string(data.y1);
CustomYLabels(end) = string(data.y2);

hm.XDisplayLabels = CustomXLabels;
hm.YDisplayLabels = CustomYLabels;

% Group Voronoi Area (product of bat areas) Heatmap (bat 1)
Z = zeros(N_points);
subplot(3,3,2);
for i=1:N_points
    for k=1:N_points
        P = cat(1, [X(i) bat_x_pos(t, 2:end)], [Y(k) bat_y_pos(t, 2:end)]);
        [A, b, V] = voronoiPolyhedronsCCW(P, [-2.8 -2.8], [2.8 2.8]);
        
        geo_mean = 1;
        for j=1:N
            V_x = V{j}(1,:);
            V_y = V{j}(2,:);
            geo_mean = geo_mean * polyarea(V_x, V_y);
        end

        Z(i,k) = geo_mean^(1/N);
    end
end

hm = heatmap(X,Y,Z);
title('Geometric Mean of all Bat Voronoi Areas');
hm.YDisplayData=flip(hm.YDisplayData);
hm.XDisplayLabels = CustomXLabels;
hm.YDisplayLabels = CustomYLabels;

% Median Group Voronoi Area (product of bat areas) Heatmap (bat 1)
subplot(3,3,3);
Z = zeros(N_points);
for i=1:N_points
    for k=1:N_points
        P = cat(1, [X(i) bat_x_pos(t, 2:end)], [Y(k) bat_y_pos(t, 2:end)]);
        [A, b, V] = voronoiPolyhedronsCCW(P, [-2.8 -2.8], [2.8 2.8]);
        
        voro_areas = [];
        for j=1:N
            V_x = V{j}(1,:);
            V_y = V{j}(2,:);
            voro_areas = [voro_areas polyarea(V_x, V_y)];
        end

        Z(i,k) = median(voro_areas);
    end
end
hm = heatmap(X,Y,Z);
title('Median of all Bat Voronoi Areas');
hm.YDisplayData=flip(hm.YDisplayData);
hm.XDisplayLabels = CustomXLabels;
hm.YDisplayLabels = CustomYLabels;

% Mean Group Voronoi Area (product of bat areas) Heatmap (bat 1)
subplot(3,3,4);
Z = zeros(N_points);
for i=1:N_points
    for k=1:N_points
        P = cat(1, [X(i) bat_x_pos(t, 2:end)], [Y(k) bat_y_pos(t, 2:end)]);
        [A, b, V] = voronoiPolyhedronsCCW(P, [-2.8 -2.8], [2.8 2.8]);
        
        voro_areas = [];
        for j=1:N
            V_x = V{j}(1,:);
            V_y = V{j}(2,:);
            voro_areas = [voro_areas polyarea(V_x, V_y)];
        end

        Z(i,k) = mean(voro_areas);
    end
end
hm = heatmap(X,Y,Z);
title('Mean of all Bat Voronoi Areas');
hm.YDisplayData=flip(hm.YDisplayData);
hm.XDisplayLabels = CustomXLabels;
hm.YDisplayLabels = CustomYLabels;

% NN Distance Heatmap (bat 1)
subplot(3,3,6);
Z = zeros(N_points);
for i=1:N_points
    for k=1:N_points
        min_dist = min(sqrt((bat_x_pos(t, 2:end)-X(i)).^2 + (bat_y_pos(t, 2:end) - Y(k)).^2));
        
        Z(i,k) = min_dist;
    end
end

hm = heatmap(X,Y,Z);
title('Nearest Neighbor Distance from Bat 1');
hm.YDisplayData=flip(hm.YDisplayData);
hm.XDisplayLabels = CustomXLabels;
hm.YDisplayLabels = CustomYLabels;

% Geometric Mean Distance Heatmap (bat 1)
subplot(3,3,7);
Z = zeros(N_points);
for i=1:N_points
    for k=1:N_points
        geomean_dist = geomean(sqrt((bat_x_pos(t, 2:end)-X(i)).^2 + (bat_y_pos(t, 2:end) - Y(k)).^2));
        
        Z(i,k) = geomean_dist;
    end
end

hm = heatmap(X,Y,Z);
title('Geometric Mean Distance from Bat 1');
hm.YDisplayData=flip(hm.YDisplayData);
hm.XDisplayLabels = CustomXLabels;
hm.YDisplayLabels = CustomYLabels;

% Mean Distance Heatmap (bat 1)
subplot(3,3,8);
Z = zeros(N_points);
for i=1:N_points
    for k=1:N_points
        mean_dist = mean(sqrt((bat_x_pos(t, 2:end)-X(i)).^2 + (bat_y_pos(t, 2:end) - Y(k)).^2));
        
        Z(i,k) = mean_dist;
    end
end

hm = heatmap(X,Y,Z);
title('Mean distance from Bat 1');
hm.YDisplayData=flip(hm.YDisplayData);
hm.XDisplayLabels = CustomXLabels;
hm.YDisplayLabels = CustomYLabels;


% Median Distance Heatmap (bat 1)
subplot(3,3,9);
Z = zeros(N_points);
for i=1:N_points
    for k=1:N_points
        median_dist = median(sqrt((bat_x_pos(t, 2:end)-X(i)).^2 + (bat_y_pos(t, 2:end) - Y(k)).^2));
        
        Z(i,k) = median_dist;
    end
end

hm = heatmap(X,Y,Z);
title('Median Distance from Bat 1');
hm.YDisplayData=flip(hm.YDisplayData);
hm.XDisplayLabels = CustomXLabels;
hm.YDisplayLabels = CustomYLabels;
end