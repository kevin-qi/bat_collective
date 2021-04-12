function [super_set] = pairwise_mirror(points)
%PAIRWISE_MIRROR Mirror all points across interbat bisection planes
%   Utility function for computing symmetrical voronoi cells
%   P := (D, N) array of bat positions. D = Number of dimensions, N =
%   number of bats
import voronoi_utils.*

N = size(points,2);
figure;
axis square;
scatter(points(1,:), points(2,:));
xlim([-3, 3]);
ylim([-3, 3]);
axis square;
hold on;

super_set = points(:,:);
for i=1:N
    for j=i+1:N
        x1 = points(1,i); y1 = points(2,i);
        x2 = points(1,j); y2 = points(2,j);
        
        del_x = x2-x1;
        del_y = y2-y1;
        
        mid_p = [[(x1+x2)/2]; [(y1+y2)/2]];
        
        p_21 = [[mid_p(1) - 1.5*del_x]; [mid_p(2) - 1.5*del_y]];
        p_12 = [[mid_p(1) + 1.5*del_x]; [mid_p(2) + 1.5*del_y]];
        
        assert(isequal((p_12 + p_21)/2, mid_p), "Mirror error, midpoints don't match!");
        
        scatter(p_21(1), p_21(2));
        scatter(p_12(1), p_12(2));
        
        super_set = [super_set p_12];
        super_set = [super_set p_21];
    end
end
        

end

