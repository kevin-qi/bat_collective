function [flight_data] = characterize_flights(data, forage_type, forage_threshold, inter_bat_threshold, verbose)
%CHARACTERIZE_FLIGHTS Characterize flight properties
% data := data of a single flight session
% forage_type := 'bowl' or 'feeders'
% forage_threshold := Euclidean distance to feeder/bowl to be considered foraging
% inter_bat_threshold := Euclidean distance to nearest bat to be considered not isolated
% verbose := bool flag to display extra figures and outputs (for debugging)


N = size(data.bflying,2);
fig_counter = 1; % Tracks number of figures to ensure new figure window

% Plot flight times of all bats 
if(verbose)
    figure(fig_counter);
    fig_counter = fig_counter + 1;
    
    for i=1:N
        subplot(N,1,i);
        plot(data.bflying(:,i));
        ylabel(data.bat_nms(i,:));
    end
    sgtitle('Flights');
end

% Minimum length of 5 bat position recordings
min_length = 999999999;
for b_num = 1:N
    if length(data.tag_data_filt{b_num}) < min_length
        min_length = length(data.tag_data_filt{b_num});
    end
end
if length(data.bflying) < min_length
    min_length = length(data.bflying)-1;
end
%disp(min_length);
% 1 at the start of flight, -1 at end of flight, 0 other wise
flight_start_end = diff(data.bflying, 1, 1); 
flight_start_end = flight_start_end(1:min_length, :);

% Get flight start end indices
for b_num=1:N
    flight_start{b_num} = find(flight_start_end(:,b_num) > 0)+1;
    flight_end{b_num} = find(flight_start_end(:,b_num) < 0);
    
    % Assert num of flights is counted correctly...
    %assert(length(flight_start{b_num}) == data.num_flights(b_num), ...
    %       'Num flights in flight_start do not match num flights in data.num_flights');
       
    %assert(length(flight_end{b_num}) == data.num_flights(b_num), ...
    %       'Num flights in flight_end do not match num flights in data.num_flights');
    
    % Assert flight start end times are in chronological order (flight
    % cannot end before it starts!)
    assert(prod(flight_end{b_num} > flight_start{b_num}), ...
           'Flight ended before flight began, check flight start/end calculations');
end

flight_data = {};
forage_count_from = zeros(N,1);
forage_count_to = zeros(N,1);
for b_num=1:N % For each bat
    flight = struct;
    for i=1:length(flight_start{b_num}) % For each flight by bat b_num
        % Flight start/end times
        flight(i).start = flight_start{b_num}(i);
        flight(i).end = flight_end{b_num}(i);     

        bflying_snippet = data.bflying(flight(i).start:flight(i).end,:);
        
        % By construction, bat b_num should be flying between flight start
        % and end.
        assert(prod(bflying_snippet(:,b_num) == 1) == 1, 'Bat b_num is not flying during supposed flight');
        
        % If flight concurrent with other flights (concurrent = 1) or 
        % solo (concurrent = 0)
        is_concurrent = max(sum(bflying_snippet, 2) > 1) == 1;
        flight(i).concurrent = is_concurrent;
        
        % Flight start (from edge or from bowl/feeder
        % TODO: Improve forage detection
        if(forage_type == 'bowl')
            r_bowl = data.bowl(1:2);
            r_bat = data.tag_data_filt{b_num}(flight(i).start, 3:4);

            dist = norm(r_bat - r_bowl);
            if(dist < forage_threshold)
                flight(i).from_edge = false;
                flight(i).from_forage = true;
                forage_count_from(b_num) = forage_count_from(b_num) + 1;
            else
                flight(i).from_edge = true;
                flight(i).from_forage = false;
            end
        end
        
        % To self flight (start end same spot)
        dist = norm(data.tag_data_filt{b_num}(flight(i).end, 3:4) - data.tag_data_filt{b_num}(flight(i).start, 3:4));
        if(dist < 0.5)
            flight(i).self = true;
        else
            flight(i).self = false;
        end
        
        % Flight end (to bat / to bowl/feeder / to solo)
        % TODO: Improve forage detection
        if(forage_type == 'bowl')
            r_bowl = data.bowl(1:2);
            r_bat = data.tag_data_filt{b_num}(flight(i).end, 3:4);
            
            dist = norm(r_bat - r_bowl);
            if(dist < forage_threshold) % To forage
                flight(i).to_edge = false;
                flight(i).to_bat = false;
                flight(i).to_forage = true;
                forage_count_to(b_num) = forage_count_to(b_num) + 1;
            else
                % nearest bat distance
                nn_dist = 99;
                for j = 1:N
                    if(j ~= b_num)
                        inter_bat_dist = norm(data.tag_data_filt{j}(flight(i).end, 3:4) - r_bat);
                        if(inter_bat_dist < nn_dist)
                            nn_dist = inter_bat_dist;
                        end
                    end
                end 
                if(nn_dist < inter_bat_threshold) % To bat
                    flight(i).to_bat = true;
                    flight(i).to_edge = false;
                    flight(i).to_forage = false;
                elseif (r_bat(1)<data.x1+0.3 || ... % To edge
                         r_bat(1)>data.x2-0.3 || ...
                         r_bat(2)<data.y1+0.3 || ...
                         r_bat(2)>data.y2-0.3)
                    flight(i).to_bat = false;
                    flight(i).to_edge = true;
                    flight(i).to_forage = false;
                end
            end
        end
        
        
    end
    flight_data{b_num} = flight; 
end

if(verbose)
    disp('Forage Counts:');
    disp(forage_count_from);
    disp(forage_count_to);
    
    red = [1, 0, 0];
    blue = [0, 0, 1];
    

    % Isolating Flights
    figure(fig_counter);
    fig_counter = fig_counter + 1;
    sgtitle('Isolating Flights');

    for b_num=1:N
        subplot(1,N,b_num);
        xlim([data.x1, data.x2]);
        ylim([data.y1, data.y2]);
        axis square;
        hold on;
        for flight_num = 1:length(flight_data{b_num})
            flight = flight_data{b_num}(flight_num);
            if(flight.to_edge)
                bat_pos = data.tag_data_filt{b_num}(flight.start:flight.end, 3:5);
                len = length(bat_pos);
                %colors_p = [linspace(red(1),blue(1),len)', linspace(red(2),blue(2),len)', linspace(red(3),blue(3),len)'];
                p = scatter(bat_pos(1,1), bat_pos(1,2), 'red');
                p = scatter(bat_pos(end,1), bat_pos(end,2), 'blue');
                % modified jet-colormap
                %cd = [uint8(colors_p*255) uint8(ones(len,1))].';
                
                %drawnow
                %set(p.Edge, 'ColorBinding', 'interpolated', 'ColorData', cd)
            end
        end
        hold off;
    end
    shg;
    
    % To Forage Flights
    figure(fig_counter);
    fig_counter = fig_counter + 1;
    sgtitle('To Forage Flights From Edge/Bat');

    for b_num=1:N
        subplot(1,N,b_num);
        xlim([data.x1, data.x2]);
        ylim([data.y1, data.y2]);
        axis square;
        hold on;
        for flight_num = 1:length(flight_data{b_num})           
            flight = flight_data{b_num}(flight_num);
            if(flight.to_forage & ~flight.from_forage)
                bat_pos = data.tag_data_filt{b_num}(flight.start:flight.end, 3:5);
                len = length(bat_pos);
                %colors_p = [linspace(red(1),blue(1),len)', linspace(red(2),blue(2),len)', linspace(red(3),blue(3),len)'];
                p = scatter(bat_pos(1,1), bat_pos(1,2), 'red');
                p = scatter(bat_pos(end,1), bat_pos(end,2), 'blue');
                % modified jet-colormap
                %cd = [uint8(colors_p*255) uint8(ones(len,1))].';
                
                %drawnow
                %set(p.Edge, 'ColorBinding','interpolated', 'ColorData', cd)
            end
        end
        hold off;
    end
    shg;
    
    % From Forage Flights
    figure(fig_counter);
    fig_counter = fig_counter + 1;
    sgtitle('From Forage Flights');

    for b_num=1:N
        subplot(1,N,b_num);
        xlim([data.x1, data.x2]);
        ylim([data.y1, data.y2]);
        axis square;
        hold on;
        for flight_num = 1:length(flight_data{b_num})           
            flight = flight_data{b_num}(flight_num);
            if(flight.from_forage)
                bat_pos = data.tag_data_filt{b_num}(flight.start:flight.end, 3:5);
                len = length(bat_pos);
                %colors_p = [linspace(red(1),blue(1),len)', linspace(red(2),blue(2),len)', linspace(red(3),blue(3),len)'];
                p = scatter(bat_pos(1,1), bat_pos(1,2), 'red');
                p = scatter(bat_pos(end,1), bat_pos(end,2), 'blue');
                % modified jet-colormap
                %cd = [uint8(colors_p*255) uint8(ones(len,1))].';
                
                %drawnow
                %set(p.Edge, 'ColorBinding','interpolated', 'ColorData', cd)
            end
        end
        hold off;
    end
    shg;
end


end

