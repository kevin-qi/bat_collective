%% Some exploratory analysis for symmetric voronoi cells

close all;

import voronoi_utils.*;

%% Equidistant points on a circle

R = 2;
theta = linspace(0,pi,36);
X = R*cos(theta);
X_mirror = R*cos(theta+pi);
Y = R*sin(theta);
Y_mirror =R*sin(theta+pi);

stat_x = X(18);
stat_x_mirr = X_mirror(18);
stat_y = Y(18);
stat_y_mirr = Y_mirror(18);

figure;
for i=1:16
    subplot(4,4,i);
    P_x = cat(2, X(i), X_mirror(i), stat_x, stat_x_mirr, 0);
    P_y = cat(2, Y(i), Y_mirror(i), stat_y, stat_y_mirr, 0);
    P = cat(1, P_x, P_y);
    disp(P);
    [A, b, V] = voronoiPolyhedronsCCW(P, [-2.8 -2.8], [2.8 2.8]);
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
        hold on
        scatter(P(1,j), P(2,j));


        text(cent_x, cent_y, string(area(j,i)));
        xlim([-3, 3]);
        ylim([-3 3]);
        axis square;            
    end
    hold off;

end



