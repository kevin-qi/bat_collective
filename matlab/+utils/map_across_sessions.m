function [result] = map_across_sessions(session_data,inputArg2)
%Map Across Sessions: Map function func across sessions
%   session_data: cell array where session_data{i} = the extracted data
%   from session i
%   func: function to be mapped

for i = 1:length(session_data) % i = index of session
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


end

