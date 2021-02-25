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
tic
%Load data, keep non-zero entries and sort according to ascending network times
file_name = dir(fullfile(cd, '*_cdp_*'));   file_name = file_name.name;
RTLS_data = load(file_name);
RTLS_data = RTLS_data(RTLS_data(:,2)~=0,:);
RTLS_data = sortrows(RTLS_data,3);
toc
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
toc