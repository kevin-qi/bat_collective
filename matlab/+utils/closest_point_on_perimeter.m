function [dist, argmin] = closest_point_on_perimeter(position, bounds, N_points)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

X = linspace(bounds.x1, bounds.x2, N_points);
Y = linspace(bounds.y1, bounds.y2, N_points);

perimeter_x = [];
perimeter_y = [];
for x=X
    perimeter_x = [perimeter_x x];
    perimeter_y = [perimeter_y bounds.y2];
end
for y=flip(Y)
    perimeter_x = [perimeter_x bounds.x2];
    perimeter_y = [perimeter_y y];
end
for x=flip(X)
    perimeter_x = [perimeter_x x];
    perimeter_y = [perimeter_y bounds.y1];
end
for y=Y
    perimeter_x = [perimeter_x bounds.x1];
    perimeter_y = [perimeter_y y];
end
del_x = position(1) - perimeter_x;
del_y = position(2) - perimeter_y;
dist = sqrt(del_x.^2 + del_y.^2);
[min_dist, argmin] = min(dist);
dist = min_dist;
end

