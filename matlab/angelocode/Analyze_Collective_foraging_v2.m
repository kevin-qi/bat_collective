%% Script for the analysis of collective foraging (v1)
%Load data
extracted_CDPfile = dir(fullfile(cd, '*extracted_*'));          load(extracted_CDPfile.name);
batdate = extracted_CDPfile.name(11:16);

%Parameters
n_tags = CDPmtdata.tags;
x2=2.8; x1=-2.8; y2=2.8;  y1=-2.8;  z1=0; z2=2.30;             %Flight volume coordinates
edges_d = {x1:(x2-x1)/10:x2 y1:(y2-y1)/10:y2};                 %Edges for density histogram
bowl = [-0.29, 0.05, 0.45];                                      %x,y,z bowl 1
Fs = 100;                                                      %resampling frequency (Hz) for common time
bat_nms = ['Dai'; 'Den'; 'Dia'; 'Dor'; 'Dum'];
bat_pairs = nchoosek(1:n_tags,2);   bat_pair_nms = [bat_nms(bat_pairs(:,1),:), '-'.*ones(length(bat_pairs),1), bat_nms(bat_pairs(:,2),:)];
bat_clr = lines(n_tags);
v_th = 0.5;                                                    %Velocity threshold (m/s)
TTL_time_diff = [21; 13; 8; 5; 4];

%Custom graded colormap
for i = 1:n_tags
    for j = 1:3
        custom_map(:,j,i) = linspace(1,bat_clr(i,j))';
    end
end

%Options
show_fig = 1;
save_data = 1;
use_r_corr = 1;
compare_c3d = 0;
savemovie = 0;

%Create analysis folder for storing the results
if save_data
    analysis_directory=fullfile(pwd,['Analysis_',batdate]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%% Calculate evenly sampled kinematic variables (r,v,a)

t = [CDPmtdata.TTL_times(1):1/Fs:CDPmtdata.TTL_times(end)]';    T = length(t);    %time vector from 1st to last TTL
r = zeros(T,3,n_tags);                                                            %3D position vector (sample,dim,id)

%Interpolate position at query time points by splining filtered data
for i = 1:n_tags
    r(:,:,i) =  csaps(tag_data_filt{1,i}(:,8), tag_data_filt{1,i}(:,[3:5])', 0.9, t)';
end

%Interpolate acceleration at query time points by splining acceleration data
if exist('tag_ac_data')
    a = zeros(T,3,n_tags);
    for i = 1:n_tags
        a(:,:,i) =  csaps(tag_ac_data{1,i}(:,8), tag_ac_data{1,i}(:,[3:5])', 1, t)';
    end
    a_abs = squeeze(vecnorm(a,2,2));    %Modulus
    a_flt = bandpass(a_abs,[7 9],100);  %Filtered at the wing-beat frequency
end

%Calculate velocity and 2D-direction of motion
v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,n_tags),v);   v_abs = squeeze(vecnorm(v,2,2));
angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));

%Detect flight epochs based on wing-beat signal
if exist('tag_ac_data')
    for i = 1:n_tags
        [up,lo] = envelope(a_flt(:,i),8,'peak');    %Envelope of the acceleration signal
        env = normalize(up - lo,'range');           %Amplitude of the envelope
        env_th = otsuthresh(histcounts(env));       %Threshold (based on Otsu method). Can be set at 0.35
        wBeats(:,i) = movsum(env>env_th,2*Fs)>Fs/3; %Euristic criterion for flight detection
        %ax(i) = subplot(n_tags,1,i);  area(t,wBeats(:,i)*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
        %plot(t,v_abs(:,i));     plot(t,r(:,1,i),t,r(:,2,i));    hold off;
    end
    %linkaxes(ax,'x');
end

%% Correct position based on wing-beats or velocity threshold
%Initialize vectors
tag_data_stat = tag_data;
r_stat = zeros(T,3,n_tags);

if exist('tag_ac_data')
    %Aggressive median filtering on original data
    for i = 1:n_tags
        tag_data_stat{1,i}(:,[3:5]) = smoothdata(tag_data_filt{1,i}(:,[3:5]),1,'movmedian',Fs*5);
        r_stat(:,:,i) =  csaps(tag_data_stat{1,i}(:,8), tag_data_stat{1,i}(:,[3:5])', 1, t)';
    end
    %Substitute median filtered data when bat is not flying
    stat_periods = repmat(~wBeats,1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
    r_corr = r_stat.*stat_periods+r.*(~stat_periods);
    r_corr = smoothdata(r_corr,1,'loess',Fs*1);
    
else
    %Aggressive median filtering on original data
    for i = 1:n_tags
        tag_data_stat{1,i}(:,[3:5]) = smoothdata(tag_data_filt{1,i}(:,[3:5]),1,'movmedian',Fs*3);
        r_stat(:,:,i) =  csaps(tag_data_stat{1,i}(:,8), tag_data_stat{1,i}(:,[3:5])', 1, t)';
    end
    %Calculate xy-velocity and smooth
    v_filt = smoothdata(squeeze( vecnorm( v(:,1:2,:),2,2)), 1, 'movmedian', Fs*3);
    %Substitute median filtered data when velocity is less than threshold
    stat_periods = repmat((v_filt < v_th),1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
    r_corr = r_stat.*stat_periods+r.*(~stat_periods);
    r_corr = smoothdata(r_corr,1,'loess',Fs*1);
end

%% Plot Position data (raw and corrected)
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,3),'.');              hold on;     sgtitle('x');
    plot(t,r_corr(:,1,i),'.');    hold off;    ylim([-3 3]);      legend('raw','corrected');
end
linkaxes(ax,'x');    xlabel('Time(s)');

figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,4),'.');              hold on;     sgtitle('y');
    plot(t,r_corr(:,2,i),'.');    hold off;    ylim([-3 3]);      legend('raw','corrected');
end
linkaxes(ax,'x');    xlabel('Time(s)');

figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,5),'.');              hold on;     sgtitle('z');
    plot(t,r_corr(:,3,i),'.');    hold off;    ylim([0 2.5]);     legend('raw','corrected');
end
linkaxes(ax,'x');    xlabel('Time(s)');

%% Correct position, velocity and angle
if use_r_corr
    r = r_corr; v_old = v_abs;
    v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,n_tags),v);   v_abs = squeeze(vecnorm(v(:,1:2,:),2,2));     %This is 2D velocity!!
    angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));
end

%% Plot velocity (raw and corrected)
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(t,v_old(:,i),'.');    hold on;     sgtitle('v');
    plot(t,v_abs(:,i),'.');    hold off;    legend('raw','corrected');
end
linkaxes(ax,'x');    xlabel('Time(s)');

%% Flight segmentation
%Detect flight starts and stops by using risetime and falltime
bflying = zeros(T,n_tags);
allsums = normalize(movsum(v_abs > v_th,Fs,1),'range');
flight_times = {};
for i = 1:n_tags
    [~,rLT,~,~,~] = risetime(allsums(:,i),'PercentReferenceLevels',[10 95]);    flight_times{i,1} = round(rLT);
    [~,fLT,~,~,~] = falltime(allsums(:,i),'PercentReferenceLevels',[10 95]);    flight_times{i,2} = round(fLT);
    
    %Ensure same number of starts/stops
    if ~isempty(flight_times{i,1})
        flight_times{i,1} = flight_times{i,1}(flight_times{i,1}<flight_times{i,2}(end,1));
        flight_times{i,2} = flight_times{i,2}(flight_times{i,2}>flight_times{i,1}(1,1));
        num_flights(i) = size(flight_times{i,2},1);
    else
        flight_times{i,1} = [];
        flight_times{i,2} = [];
        num_flights(i) = 0;
    end
    
    for f = 1:num_flights(i)
        bflying(flight_times{i,1}(f,1):flight_times{i,2}(f,1),i) = 1;
    end
end

%Define staionary periods when the bat is not flying
stat_periods = repmat(~bflying,1,1,3);   stat_periods = permute(stat_periods,[1 3 2]);
r_qt = r;   r_qt(~stat_periods)=nan;     angle(squeeze(stat_periods(:,1,:))) = nan;

%% Feeding attempts
num_feeds = zeros(n_tags,1);
bfeeds = zeros(T,n_tags);
for i = 1:n_tags
    feed_times{i,1} = [];
    feed_times{i,2} = [];
    
    for f = 1:num_flights(i)
        if vecnorm(squeeze(r(flight_times{i,2}(f,1),1:2,i))-bowl(1:2),2,2) < 0.3
            num_feeds(i,1) = num_feeds(i,1)+1;
            feed_times{i,1}(num_feeds(i,1),1) = flight_times{i,2}(f,1);
            feed_times{i,2}(num_feeds(i,1),1) = flight_times{i,1}(f+1,1);
        end
    end
    
    for f = 1:num_feeds(i)
        bfeeds(feed_times{i,1}(f,1):feed_times{i,2}(f,1),i) = 1;
    end
end

%% Social Network Graph
%Adjacency matrix for proximity
A = zeros(n_tags);
for i = 1:length(bat_pairs)
    inter_bat = [];     inter_bat = vecnorm(r_qt(:,:,bat_pairs(i,1))-r_qt(:,:,bat_pairs(i,2)),2,2);     inter_bat = inter_bat(~isnan(inter_bat));
    A(bat_pairs(i,1),bat_pairs(i,2)) = nnz(inter_bat<0.3)/T;
end

%Create Graph
G = graph(A,cellstr(bat_nms),'upper');

%% Bat landing on a bat
for i = 1:n_tags
    if ~isempty(flight_times{i,1})
        for n =flight_times{i,2}
            land_table{i,1} = squeeze(vecnorm(r(n,1:2,i)-r(n,1:2,:),2,2)<0.15);
        end
        land_table{i,1}(:,i) = 0;
    else
        land_table{i,1} = [];
    end
end

%% Probabilities and angles when flying together

p_fly1 = sum(bflying,1)./T;
heading_diff = nan(T,length(bat_pairs));
for i = 1:length(bat_pairs)
    %Probabilities
    p_fly2(i,1) = sum(bflying(:,bat_pairs(i,1)).*bflying(:,bat_pairs(i,2)),1)/T;
    p_fly2(i,2) = p_fly1(bat_pairs(i,1))*p_fly1(bat_pairs(i,2));
    
    %Angles
    cpled = find(bflying(:,bat_pairs(i,1)).*bflying(:,bat_pairs(i,2)));
    heading_diff(cpled,i) = rad2deg(angdiff(angle(cpled,bat_pairs(i,1)),angle(cpled,bat_pairs(i,2))));
end

%% Effect of time-shifts
lags = [round(-4*Fs):1:round(4*Fs)];
p_fly_sh = zeros(length(lags),length(bat_pairs));
for i = 1:length(bat_pairs)
    m = 1;
    for k = lags
        bflying_sh = circshift(bflying(:,bat_pairs(i,1)),k);
        p_fly_sh(m,i) = sum(bflying_sh.*bflying(:,bat_pairs(i,2)),1)/T;
        m = m+1;
    end
end

%% Make a few figures
if show_fig
    
    %FIG: Scatter plot all bats
    figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.5]);
    for i=1:n_tags
        subplot(131);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim([x1 x2]); ylim([y1 y2]); zlim([z1 z2]);  title('3D view');                 hold on;
        subplot(132);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim([x1 x2]); ylim([y1 y2]); zlim([z1 z2]);  title('Top view');   view(0,90);  hold on;
        subplot(133);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim([x1 x2]); ylim([y1 y2]); zlim([z1 z2]);  title('Door view');  view(-30,0); hold on;
    end
    hold off;
    
    %FIG: Trajectories individual bats
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        subplot(1,5,i);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'-','Color', bat_clr(i,:));  xlim([x1 x2]); ylim([y1 y2]); zlim([z1 z2]);  title(bat_nms(i,:));  view(0,90);
    end
    
    %FIG: Density Histograms 3D
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        subplot(1,5,i); hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','FaceColor',bat_clr(i,:));  xlim([x1 x2]); ylim([y1 y2]);  title(bat_nms(i,:));
    end
    
    %FIG: Density Histograms heat-map
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        sbpt = subplot(1,5,i);
        hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','edgecolor','none','FaceColor','interp');
        xlabel('x');
        xlim([x1 x2]); ylim([y1 y2]);   title(bat_nms(i,:));  view(0,90);  colormap(sbpt,custom_map(:,:,i)); % Change color scheme
    end        
    
    %FIG: Velocity and flight segmentation
    figure();       set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:n_tags
        ax(i) = subplot(n_tags,1,i);
        area(t,bflying(:,i)*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
        area(t,wBeats(:,i)*-1,'FaceAlpha',0.3,'LineStyle','none');
        plot(t,v_abs(:,i),'.','Color', bat_clr(i,:));     plot(t,r(:,1,i),'k--');  ylabel('Velocity (m/s)');     hold off;
        legend('Fly','Wing-B','Vel','x(m)');
        title([num2str(num_flights(i)) ' flights']);
    end
    linkaxes(ax,'x');   xlabel('Time(s)');
    
    %FIG: Few example flights and direction of motion
    for i = 1:n_tags
        figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        nb = 1; sgtitle(bat_nms(i,:));
        for n = randperm(min(num_flights),min([15,num_flights]))
            chunk = flight_times{i,1}(n,1):flight_times{i,2}(n,1);
            subplot(3,5,nb);     plot(r(chunk,1,i),r(chunk,2,i),'k','LineWidth',3);    hold on;
            quiver(downsample(r(chunk,1,i),20),...
                downsample(r(chunk,2,i),20),...
                downsample(v_abs(chunk,i),20).*downsample(cos(angle(chunk,i)),20),...
                downsample(v_abs(chunk,i),20).*downsample(sin(angle(chunk,i)),20),1.5,'k');
            hold off;
            xlim([x1 x2]); ylim([y1 y2]);
            nb = nb+1;
        end
    end
    
%     %FIG: Correlation velocities
%     figure();   corr_v = corr(v_filt);    imagesc(corr_v); title('Velocity (filtered) correlations');
%     xticks([1:n_tags]); xticklabels(cellstr(bat_nms));
%     yticks([1:n_tags]); yticklabels(cellstr(bat_nms));
    
    %FIG: Inter-bat distances (when both bats are stationary)
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(bat_pairs)
        subplot(2,5,i); histogram(vecnorm(r_qt(:,:,bat_pairs(i,1))-r_qt(:,:,bat_pairs(i,2)),2,2),'edgecolor','none');
        xlabel('Inter-bat distance (m)'); title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:)]);
        yticks([]); set(gca,'XScale','log')
    end
    
    %FIG: Feeding attemps
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    for i=1:n_tags
        cx(i) = subplot(n_tags,1,i);
        if ~isempty(feed_times{i,1})
            area(t,bfeeds(:,i),'FaceColor', bat_clr(i,:));  ylabel('Bat at the bowl');
        end
        title([num2str(num_feeds(i)) ' feeds']);
    end
    linkaxes(cx,'x');    xlabel('Time(s)');
    
    %FIG: Network-graph
    %Plot the edges such that their line width is proportional to their weight
    figure();
    G.Edges.LWidths = 20*G.Edges.Weight; %7*G.Edges.Weight/max(G.Edges.Weight);
    p_NG = plot(G);
    p_NG.LineWidth = G.Edges.LWidths;  p_NG.MarkerSize = 10;  p_NG.NodeLabelColor = bat_clr;
    p_NG.NodeFontSize = 15;    p_NG.NodeFontWeight = 'bold';
    
    %FIG: Trajectories individual bats on top of histogram of the others
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        colormap('copper') % Change color scheme
        subplot(1,5,i);
        hist3(reshape(permute(r(:,1:2,find(1:n_tags ~= i)),[1 3 2]),[],2),'edges',edges_d,'CdataMode','auto','edgecolor','none','FaceColor','interp');
        xlabel('x');    view(0,90);   hold on;
        plot3(r(:,1,i),r(:,2,i),1e9.*r(:,3,i),'-','Color', 'w');  xlim([x1 x2]); ylim([y1 y2]);   title(bat_nms(i,:));  view(0,90);
        hold off;
    end
    
    %FIG: Probabilities of flying together
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);
    scatter(p_fly2(:,2),p_fly2(:,1),'filled');
    xlabel('P coupled fly (independent)');  ylabel('P coupled fly (real)');
    text(p_fly2(:,2),p_fly2(:,1),bat_pair_nms,'VerticalAlignment','bottom','HorizontalAlignment','right')
    axis square;    xlim(ylim);     hold on;    h = refline(1,0);   hold off;   h.Color = 'k';  h.LineStyle = '--';
    
    %FIG: Effect of time shifts
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(bat_pairs)
        subplot(2,5,i); plot(lags./Fs,p_fly_sh(:,i));   hold on;    xline(0);   hold off;
        xlabel('Time shift 1st bat (s)'); title([bat_nms(bat_pairs(i,1),:) '-' bat_nms(bat_pairs(i,2),:)]);
        ylabel('Probability of pair flight');
        ylim([0 0.04]);
    end
    
    %FIG: Bat landing on a bat
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i =1:n_tags
        if ~isempty(flight_times{i,1})
            subplot(1,n_tags,i);
            bar(sum(land_table{i,1}),'facecolor',bat_clr(i,:));
        end
    end
    
    %FIG: Percentages of time flying together
    figure();   
    histogram(sum(bflying,2),'Normalization','probability');
    ylabel('Probability');  xlabel('Bats simult. flying');
    
end

%% Save figures
if save_data
    figHandles = findall(0,'Type','figure');
    for i = 1:numel(figHandles)
        saveas(figHandles(i),[analysis_directory, '/', batdate '_figure' num2str(i) '.png']);
    end
    close all;
end

%% Save data
if save_data
    Table=NaN(n_tags,3);
    Table(:,1)= (1:1:n_tags)';
    Table(:,2)= (num_flights)';
    Table(:,3)= (num_feeds)';
    dlmwrite('Table_j.txt',Table,'delimiter','\t','precision',4)
    save([analysis_directory,'/Analysis_', batdate, '.mat']);
end

%% Save movie with bat trajectories (45 min)
if savemovie
    close all;
    figure();       set(gcf, 'units','normalized','outerposition',[0.3 0.3 0.21 0.4]);
    n=1;
    tic;
    %movie = struct('cdata',repmat({uint8(zeros(420,560,3))},1,100),'colormap',[]);
    for j=1:20:270000   %length(r)
        for i = 1:n_tags
            plot(r_corr(j,1,i),r_corr(j,2,i),'.','MarkerFaceColor', bat_clr(i,:),'MarkerSize',60);  hold on;
        end
        plot(b2(1),b2(2),'o','MarkerSize',10);
        hold off;
        xlim([x1 x2]); ylim([y1 y2]);   xticks([]); yticks([]);
        drawnow();
        movie(n) = getframe(gcf) ;
        n=n+1;
    end
    
    % create the video writer with 1 fps
    writerObj = VideoWriter([batdate '-myVideo.avi']);
    writerObj.FrameRate = 10;
    open(writerObj);
    for i=1:length(movie)
        writeVideo(writerObj, movie(i));
    end
    close(writerObj);
    toc
    close all;
end