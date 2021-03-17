%% Script for the analysis of collective flight behavior

data = load('data/Analysis_210222/Analysis_210222.mat');
bflying = data.bflying;
N = min(size(bflying));
fs = 100;
total_flying = sum(bflying, 2);
%plot(total_flying);
%shg;

time = linspace(1, length(total_flying)/100, length(total_flying));
cum_flights = cumsum(total_flying);

% Runtime flags (0 to skip that block)
fig2_flag = 0;
fig3_flag = 0;
fig4_flag = 0;
fig5_flag = 0;
fig6_flag = 0;
fig7_flag = 0;
fig8_flag = 0;
fig9_flag = 0;
fig10_flag = 0;
fig11_flag = 1;
fig12_flag = 1;
fig13_flag = 1;
fig14_flag = 1;
fig15_flag = 1;

%% Cumulative Average Flight Time
if(fig2_flag)
    figure(2);
    time = linspace(1, length(total_flying)/100, length(total_flying));
    plot(time, cum_flights/100);
    xlabel('time (s)');
    ylabel('cumulated flight time (s)');
    title('Cumulative Total Flight Time');
    shg;
end


%% Cumulative Flight Time per Bat
if(fig3_flag)
    figure(3);
    plot(time, cumsum(bflying(:,1))/100, ...
         time, cumsum(bflying(:,2))/100, ...
         time, cumsum(bflying(:,3))/100, ...
         time, cumsum(bflying(:,4))/100, ...
         time, cumsum(bflying(:,5))/100);
    title('Cumulative Flight Time per Bat');
    xlabel('time (s)');
    ylabel('cumulated flight time (s)');
    shg;
end


%% Power law fit for cumulative flight time
if(fig4_flag)
    figure(4);
    for i = 1:1
        cum_flight_i = cumsum(bflying(:,i))/100;
        [f, gof_info] = fit(time.',cum_flight_i,'b*x^m', 'StartPoint', [1, 1]);
        cum_flight_power(:, i) = f.b*time.^f.m;
        plot(time, cum_flight_power(:, i), time, cum_flight_i);
        xlabel('time (s)');
        ylabel('cumulated flight time (s)');
        title('Log-log Scaled Cumulative Flight Time');
        set(gca, 'YScale', 'log', 'XScale', 'log');
        hold on
    end
    hold off
    
    shg;
end




%% Power law distribution
if(fig5_flag)
    figure(5);
    for i = 1:N
        plot(diff(cum_flight_power(100:length(total_flying),i)/max(cum_flight_power(:,i))));
        hold on
    end
    hold off
    shg;
end

%% Flight Transitions
if(fig6_flag)
    figure(6);
    for i = 1:N
        plot((diff(bflying(:,i))>0)+2*i);
        ylim([-0.5,11.5]);

        hold on
    end
    hold off
    shg;
end

%% Cumulative Flight Transitions
if(fig7_flag)
    figure(7);
    for i = 1:N
        plot(time(1,1:length(time)-1), cumsum((diff(bflying(:,i))>0)));

        hold on
    end
    xlabel('time (s)');
    ylabel('# of flights');
    title('Cumulated # of Flights');
    hold off
    shg;
end

%% Inter-flight-Intervals (IFI)
if(fig8_flag)
    figure(8);
    num_bins = 20;
    bin_width = 20;
    bins = linspace(bin_width, bin_width*num_bins, num_bins);
    for i = 1:N
        subplot(1,5,i);
        
        flight_indicator = (diff(bflying(:,i))>0);
        time = linspace(1,length(flight_indicator)/100, length(flight_indicator));
        ifi = diff(time(flight_indicator)).';
        
        ifi_hist = histogram(ifi, bins, 'Normalization', 'pdf');
        
        hold on
        
        [lambda] = expfit(ifi_hist.BinCounts.'/length(ifi));
        flight_rate(i) = lambda; % IFI exponential distribution rate parameter for bat i
        
        plot(bins, lambda*exp(-lambda*(bins)), 'LineWidth', 2);
        ylim([0 0.02]);
        title('Bat ' + string(i) + ' rate = ' + string(lambda))
        hold off
    end
    
    shg;
end

%% -------------------------------------------------
%% In Silico Realizations of Poisson Arrival Process 
%% Based on IFI's calculated above
%% -------------------------------------------------

%% In Silico Flight Times
if(fig9_flag)
    figure(9);
    
    T = 3600; % 1 hr duration
    delta = 1/100; % 100hz sampling rate
    N_steps = T/delta; % Number of samples
    event = zeros(N, N_steps);
    for i = 1:N
        event(i, :) = zeros(1, N_steps);
        R = rand(size(event(i, :)));
        event(i, R<flight_rate(i)*delta) = 1; % Flight occurs if R < lambda*delta
        disp(size(event(i, :)))
        plot(event(i, :) + 2*i);
        ylim([-0.5 11.5]);
        hold on
    end
    hold off
    shg;
    
end

%% In Silico Cumulative Flight Times
if(fig10_flag)
    figure(10);
    
    for i = 1:N
        plot(cumsum(event(i, :)));
        hold on
    end   
    hold off
    shg;
end

%% Power law fit for cumulative flights
if(fig11_flag)
    figure(11);
    
    for i = 1:N
        flights_i = cumsum((diff(bflying(:,i))>0));
        time = linspace(1,length(flights_i), length(flights_i));
        [f, gof_info] = fit(time.',flights_i,'m*x^b', 'StartPoint', [1, 1]);
        flights_power(:, i) = f.m*time.^f.b;
        disp(f.m);
        disp(f.b);
        plot(time, flights_power(:, i), time, flights_i);
        xlabel('time (s)');
        ylabel('cumulated flights');
        title('Cumulated Flights Power Law Curves');
        hold on
    end
    hold off
    shg;
end

