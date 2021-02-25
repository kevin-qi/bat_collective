%% Extract tracking data from Ciholas RTLS acquisition
% Data are recorded through a python script and saved in a ASCII file
% Data columns (comma delimited) are defines as follows:

% Column 1 = pos(1), sync(2) or acc(0) type index
% Column 2 = device Serial Number
% Column 3 = Network time (each tic is ~15.65 ps)
% Column 4,5,6 = x,y,z (position in mm or acceleration-sync)
% Column 7 = signal quality for pos or scale for acc
% Column 8 = number of receiving anchors for a tag or '0' for acc/sync

% The sync signal is driven by a TTL:
% 50 ms duration happening every 21, 13, 8, 5, 4s
%-----------------------------------------------------------------
clear all;  close all;

%Load data, keep non-zero entries and sort according to ascending network times
file_name = dir(fullfile(cd, '*_cdp_*'));   file_name = file_name.name;
RTLS_data = load(file_name);
RTLS_data = RTLS_data(RTLS_data(:,2)~=0,:);
RTLS_data = sortrows(RTLS_data,3);

%Parameters
use_sync = 1;
rec_duration = 3600;                                                                            %approx rec duration (s)
ntw_time = 15.65e-12;                                                                           %network tic interval (s)
tags_SN = [17106917; 17106934; 17107055; 17106969; 17107100; 17107128];                         %tag serial numbers
sync_SN = 17040920; %17106963;                                                                  %sync tag serial number
n_tags = length(tags_SN);
TTL_time_diff = [21; 13; 8; 5; 4];                                                              %TTL delays in s
TTL_abs_times = [0; cumsum(repmat(TTL_time_diff,round(rec_duration*2/sum(TTL_time_diff)),1))];
CDPmtdata.Fs = 100;                                                                             %Acquisition frequency CDP (Hz)

%Extract data and log CDP metadata, start from sync
sync_data = RTLS_data(RTLS_data(:,1)==2 & RTLS_data(:,2)==sync_SN,2:end);
CDPmtdata.sync = ~isempty(sync_data);
CDPmtdata.sync_duration = (sync_data(end,2)-sync_data(1,2))*ntw_time/60;
CDPmtdata.sync_samples = length(sync_data);
%Tag Data
j = 0;
for i = 1:n_tags
    if ~isempty(find(RTLS_data(:,2)==tags_SN(i)))
        j = j+1;
        tag_data{j} = RTLS_data(RTLS_data(:,1)==1 & RTLS_data(:,2)==tags_SN(i),2:end);
        tag_data{1,j}(:,[3:5]) = tag_data{1,j}(:,[3:5])/1000;
        CDPmtdata.tag_duration(j) = (tag_data{1,j}(end,2)-tag_data{1,j}(1,2))*ntw_time/60;
        CDPmtdata.tag_samples(j) = length(tag_data{1,j});
        
        if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(i))))
            tag_ac_data{j} = RTLS_data(RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(i),2:end);
            indexes = [false(size(tag_ac_data{1,j}(:,1:2))), tag_ac_data{1,j}(:,[3:5])>double(intmax('int32'))];
            tag_ac_data{1,j}(indexes) = tag_ac_data{1,j}(indexes)-double(intmax('uint32'));
            tag_ac_data{1,j}(:,[3:5]) = tag_ac_data{1,j}(:,[3:5])*2/double(intmax('int32'));
            CDPmtdata.tag_ac_duration(j) = (tag_ac_data{1,j}(end,2)-tag_ac_data{1,j}(1,2))*ntw_time/60;
            CDPmtdata.tag_ac_samples(j) = length(tag_ac_data{1,j});
        end
    end
end

n_tags = length(tag_data);
CDPmtdata.tags = n_tags;
disp(CDPmtdata);

%% Interpolate time stamps (s) according to sync signal

if use_sync
    %Detect network times corresponding to peaks in the sync signal
    %[~, TTL_network_times] = findpeaks(normalize([0; diff(sync_data(:,3))],'zscore'),sync_data(:,2),'MinPeakHeight',2,'MinPeakDistance',2/ntw_time);
    [~, TTL_network_times] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,2),'MinPeakHeight',1,'MinPeakDistance',2/ntw_time);
    
    %Check if all the TTL were correctly detected
    [~, control] = findpeaks([0; diff(TTL_network_times)*ntw_time],'MinPeakHeight',20);
    if unique(diff(control))==5
        disp('All TTLs were detected, alignment OK');
    else
        disp('Undetected TTLs, check alignment process');
    end
    subplot(211);   plot(diff(control),'.');    ylabel('TTL count');     xlabel('#TTL pair separated by 21s');
    subplot(212);   findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,2),'MinPeakHeight',1,'MinPeakDistance',2/ntw_time);   xlabel('Packet#');     ylabel('Accelerometer (a.u.)');
    
    %Assign to each network time the corresponding absolute time (s),
    %starting from the first TTL
    TTL_abs_times = TTL_abs_times(1:length(TTL_network_times));
    CDPmtdata.TTL_times = TTL_abs_times;
    cdp_t = interp1(TTL_network_times,TTL_abs_times,sync_data(:,2),'linear','extrap');
    sync_data(:,8) = cdp_t;
    for j = 1:n_tags
        tag_data{1,j}(:,8) = interp1(sync_data(:,2),cdp_t,tag_data{1,j}(:,2),'linear','extrap');
        if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(j))))
            tag_ac_data{1,j}(:,8) = interp1(sync_data(:,2),cdp_t,tag_ac_data{1,j}(:,2),'linear','extrap');
        end
    end
else
    sync_data(:,8) = (sync_data(:,2)-sync_data(1,2)).*ntw_time;
    CDPmtdata.TTL_times = TTL_abs_times(TTL_abs_times<sync_data(end,8));
    for j = 1:n_tags
        tag_data{1,j}(:,8) = (tag_data{1,j}(:,2)-sync_data(1,2)).*ntw_time;
        if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(j))))
            tag_ac_data{1,j}(:,8) = (tag_ac_data{1,j}(:,2)-sync_data(1,2)).*ntw_time;
        end
    end
end

%% Filter position data
tag_data_filt = tag_data;
for i = 1:n_tags
    %tag_data_filt{1,i}(:,[3:5]) = smoothdata(tag_data{1,i}(:,[3:5]),1,'movmedian',CDPmtdata.Fs*2);
    %tag_data_filt{1,i}(:,[3:5]) = smoothdata(tag_data{1,i}(:,[3:5]),1,'gaussian',CDPmtdata.Fs*1);
    %tag_data_filt{1,i}(:,[3:5]) = lowpass(tag_data_filt{1,i}(:,[3:5]),0.1,CDPmtdata.Fs);
    tag_data_filt{1,i}(:,[3:5]) = smoothdata(tag_data_filt{1,i}(:,[3:5]),1,'loess',CDPmtdata.Fs*2);
end

%Save data
if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(1))))
    save(['extracted_', file_name(1:end-4), '.mat'],'tag_data','sync_data','CDPmtdata','tag_data_filt','tag_ac_data');
else
    save(['extracted_', file_name(1:end-4), '.mat'],'tag_data','sync_data','CDPmtdata','tag_data_filt');
end

%% Plot Raw and Filtered position data
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,3),'.');              hold on;     sgtitle('x');
    plot(tag_data_filt{1, i}(:,8), tag_data_filt{1, i}(:,3),'.');    hold off;    ylim([-3 3]);  legend('raw','filt');
end
linkaxes(ax,'x');    xlabel('Time(s)');

figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,4),'.');              hold on;     sgtitle('y');
    plot(tag_data_filt{1, i}(:,8), tag_data_filt{1, i}(:,4),'.');    hold off;    ylim([-3 3]);  legend('raw','filt');
end
linkaxes(ax,'x');    xlabel('Time(s)');

figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,5),'.');              hold on;     sgtitle('z');
    plot(tag_data_filt{1, i}(:,8), tag_data_filt{1, i}(:,5),'.');    hold off;    ylim([0 2.5]); legend('raw','filt');
end
linkaxes(ax,'x');    xlabel('Time(s)');

%% Sanity check on aquisition intervals and position values
%Sync
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   sgtitle('Sync data');
subplot(121);     plot(diff(sync_data(:,2))*ntw_time);         xlabel('Packet#');  ylabel('Network Time Difference(s)');
subplot(122);     histogram(diff(sync_data(:,2))*ntw_time);    set(gca,'YScale','log');   xlabel('Network Time Difference(s)');  ylabel('Counts');
disp(['Sync: Fs (network time): ' num2str(1/mean(diff(sync_data(:,2))*ntw_time)) 'Hz']);
disp(['Sync: Fs: ' num2str(1/mean(diff(sync_data(:,8)))) 'Hz']);

%Tags
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   sgtitle('Tags data');
for i=1:n_tags
    ax(i) = subplot(n_tags,2,2*i-1);   plot(diff(tag_data{1, i}(:,2))*ntw_time);          ylabel('Network Time Difference(s)');          xlabel('Packet#');
    bx(i) = subplot(n_tags,2,2*i);     histogram(diff(tag_data{1, i}(:,2)*ntw_time),1);   set(gca,'YScale','log');     ylabel('Counts'); xlabel('Network Time Difference(s)');
    disp(['Tag' num2str(i) ': Fs (network time): ' num2str(1/mean(diff(tag_data{1, i}(:,2))*ntw_time)) 'Hz']);
    disp(['Tag' num2str(i) ': Fs: ' num2str(1/mean(diff(tag_data{1, i}(:,8)))) 'Hz']);
end
linkaxes(ax,'x');       linkaxes(bx,'x');

%% Tags acceleration
if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(1))))
    
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   sgtitle('Tags data');
    for i=1:n_tags
        ax(i) = subplot(n_tags,2,2*i-1);   plot(diff(tag_ac_data{1, i}(:,2))*ntw_time);          ylabel('Network Time Difference(s)');          xlabel('Packet#');
        bx(i) = subplot(n_tags,2,2*i);     histogram(diff(tag_ac_data{1, i}(:,2)*ntw_time),1);   set(gca,'YScale','log');     ylabel('Counts'); xlabel('Network Time Difference(s)');
        disp(['Tag acceleration' num2str(i) ': Fs (network time): ' num2str(1/mean(diff(tag_ac_data{1, i}(:,2))*ntw_time)) 'Hz']);
        
    end
    linkaxes(ax,'x');       linkaxes(bx,'x');
    
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
    for i=1:n_tags
        ax(i) = subplot(n_tags,1,i);
        plot(tag_ac_data{1, i}(:,8), tag_ac_data{1, i}(:,3),'.');   sgtitle('x');   
        %hold on;    refline(0,double(intmax('int32')));     hold off;
    end
    linkaxes(ax,'x');    xlabel('Time(s)');
    
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
    for i=1:n_tags
        ax(i) = subplot(n_tags,1,i);
        plot(tag_ac_data{1, i}(:,8), tag_ac_data{1, i}(:,4),'.');   sgtitle('y');  
        %hold on;    refline(0,double(intmax('int32')));     hold off;
    end
    linkaxes(ax,'x');    xlabel('Time(s)');
    
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);   %x_axis
    for i=1:n_tags
        ax(i) = subplot(n_tags,1,i);
        plot(tag_ac_data{1, i}(:,8), tag_ac_data{1, i}(:,5),'.');   sgtitle('z');   
        %hold on;    refline(0,double(intmax('int32')));     hold off;
    end
    linkaxes(ax,'x');    xlabel('Time(s)');
    
end
