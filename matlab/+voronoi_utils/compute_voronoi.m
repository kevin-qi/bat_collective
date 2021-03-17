function [voronoi_areas] = compute_voronoi(data, run_tests, show_figures)
%compute_voronoi Compute voronoi areas of bats at each unique time point
%   Unique time points are defined as times with different bat positions
%   after a flight. 
%   ----------------
%   data : processed rtls data from Analyze_Collective_Foraging_*.m scripts
%   run_tests : Bool flag to run tests (monte carlo area calculations etc.)
%   show_figures : Bool flag to generate figures + movies (Takes a
%   significant amount of time)

import voronoi_utils.*;

N = 5;
T = 355000;
t_start = 10000;

% Filtered Position
tag_data_filt = data.tag_data_filt(1:5);
x_pos = zeros(N,T);
y_pos = zeros(N,T);
for i = 1:N
    x_pos(i,:,:) = tag_data_filt{1, i}(1:T, 3).';
    y_pos(i,:,:) = tag_data_filt{1, i}(1:T, 4).';
end


% Flight room bounds
bounds = [data.x1 data.y1;
          data.x2 data.y1;
          data.x1 data.y2;
          data.x2 data.y2;];

%% Get position configurations during stationary periods
%% without duplicates (1 config per stationary segment)
is_dupe = false;
unique_times = [];
t0 = 1; % start of stat period
t1 = 1; % end of stat period
for t = t_start:T
    if sum(data.bflying(t,:)) == 0 % No bats are flying
        if(~is_dupe)
            t0 = t; % Start of stat period
            is_dupe = true;
        end
    else % reset is_dupe flag and keep parsing till next stationary period
        
        if(is_dupe)
            t1 = t;
            if(t1 - t0 > 1000)
                unique_times(end + 1) = floor((t0 + t1)/2);
                is_dupe = false;
            end
        end
    end
end

if(run_tests)
    figure(4)
    plot(data.stat_periods(:,1));
    hold on;
    time_points = zeros(length(data.bflying), 1);
    time_points(unique_times) = 1;
    plot(time_points);
    ylim([0 2]);
    hold off;
    shg;


    disp(size(unique_times));
    % Monte Carlo Voronoi Area Estimation
    bat_areas = zeros(length(unique_times),5);
    for i = 1:length(unique_times)
        t = unique_times(i);
        counts = zeros(5,1);
        N_points = 1000000;
        X = rand(N_points, 1)*(data.x2 - data.x1) - data.x2;
        Y = rand(N_points, 1)*(data.y2 - data.y1) - data.y2;
        R = cat(2,X,Y);
        x = x_pos(:,t);
        y = y_pos(:,t);
        r = cat(2,x,y);

        pair_dist = pdist2(R,r);
        [M,I] = min(pair_dist,[],2);

        for j = 1:N
            counts(j) = sum(I(:) == j);
        end
        A = (data.x2 - data.x1)*(data.y2 - data.y1);
        c_fractions = counts/N_points;
        bat_areas(i,:) = A*c_fractions;
    end
    
    figure(1);
    scatter(X,Y);
    title('Monte Carlo Sampling Points');
    shg;
    figure(2);
    hist(X);
    title('Monte Carlo X position distribution');
    shg;
    figure(3);
    hist(Y);
    title('Monte Carlo Y position distribution');
    shg;
end

%% Compute Voronoi Cell Areas
for i=1:length(unique_times)
    t = unique_times(i);
    P = cat(2, x_pos(:,t),y_pos(:,t));
    [A, b, V] = voronoiPolyhedronsCCW(P.', [-2.8 -2.8], [2.8 2.8]);
    for j=1:length(V)
        % X & Y coordinates of the j-th voronoi cell vertices
        V_x = V{j}(1,:);
        V_y = V{j}(2,:);

        area(j,i) = polyarea(V_x, V_y);
    end
end

if(show_figures)
    figure(5);
    F = getframe;
    for i=1:length(unique_times(1:60))
        t = unique_times(i);
        P = cat(2, x_pos(:,t),y_pos(:,t));
        [A, b, V] = voronoiPolyhedronsCCW(P.', [-2.8 -2.8], [2.8 2.8]);
        for j=1:length(V)
            
            % X & Y coordinates of the j-th voronoi cell vertices
            V_x = V{j}(1,:);
            V_y = V{j}(2,:);
            
            area(j,i) = polyarea(V_x, V_y);
            
            % Compute centroid
            cent_x = mean(V_x);
            cent_y = mean(V_y);
            
            % Complete the polygon by appending the first entry (for
            % plotting purposes)
            V_x(end+1) = V_x(1); 
            V_y(end+1) = V_y(1);
            
            plot(V_x, V_y);
            hold on;
            scatter(P(j,1), P(j,2));
            
            
            text(cent_x, cent_y+0.1, string(area(j,i)));
            if(run_tests)
                text(cent_x, cent_y-0.1, string(bat_areas(i,j)));
            end
            ylim([data.y1-0.2, data.y2+0.2]);
            xlim([data.x1-0.2, data.x2+0.2]);
            axis square;            
        end
        F(length(F)+1) = getframe(gcf);
        hold off;
    end

    writer = VideoWriter('voronoi');
    writer.FrameRate = 1;
    open(writer);
    for frame=F(2:length(F))
        writeVideo(writer, frame);
    end
    close(writer);
end

voronoi_areas = area;
end

