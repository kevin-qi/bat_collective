function [pairwise_dist] = pairwise_distance(positions)
%pairwise_distance Calculate pairwise distances between bats

N = size(positions, 2);

for i=1:N
    for j=1:N
        r1 = positions{i};
        r2 = positions{j};
        if(i==j)
            pairwise_dist{i}(:, j) = linspace(inf, inf, length(r1));
        else
            pairwise_dist{i}(:, j) = vecnorm(r2-r1,2,2);
        end
    end
end

end

