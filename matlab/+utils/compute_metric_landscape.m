function [Z] = compute_metric_landscape(bat_pos, b_num, metric, bounds, N_points)
%COMPUTE_METRIC_LANDSCAPE Compute the landscape of b_num given the positions of
%other bats (X) based on metric with a sampling density set by N_points in
%a room set by bounds

import voronoi_utils.*;

x1 = bounds.x1;
x2 = bounds.x2;
y1 = bounds.y1;
y2 = bounds.y2;

N = length(bat_pos); % number of bats

N_points = 40;
X = linspace(x1, x2, N_points);
Y = linspace(y1, y2, N_points);
[A,B] = meshgrid(X, Y);
XY = [A(:) B(:)];

Z = zeros(N_points);
for i=1:N_points
    for k=1:N_points
        bat_pos(1, b_num) = X(i);
        bat_pos(2, b_num) = Y(k);
        [A, b, V] = voronoiPolyhedronsCCW(bat_pos, [x1 y1], [x2 y2]);
        
        other_bat_ind = linspace(1,N,N);
        other_bat_ind = other_bat_ind(other_bat_ind ~= b_num);
        bat_x_pos = bat_pos(1, other_bat_ind);
        bat_y_pos = bat_pos(2, other_bat_ind);
        
        if strcmp(metric, 'voro-geo-mean')
            geo_mean = 1;
            for j=1:N
                V_x = V{j}(1,:);
                V_y = V{j}(2,:);
                geo_mean = geo_mean * polyarea(V_x, V_y);
            end

            Z(N_points-k+1,i) = geo_mean^(1/N);
        elseif strcmp(metric, 'voro')
            V_x = V{1}(1,:);
            V_y = V{1}(2,:);

            Z(N_points-k+1,i) = polyarea(V_x, V_y);
        elseif strcmp(metric, 'nearest-neighbor')
            min_dist = min(sqrt((bat_x_pos-X(i)).^2 + (bat_y_pos - Y(k)).^2));
        
            Z(N_points-k+1,i) = min_dist;
        elseif strcmp(metric, 'mean-dist')
            mean_dist = mean(sqrt((bat_x_pos-X(i)).^2 + (bat_y_pos - Y(k)).^2));

            Z(N_points-k+1,i) = mean_dist;
        elseif strcmp(metric, 'test')
            Z(N_points-k+1,i) = i+k;
        end
            
    end
end
end

