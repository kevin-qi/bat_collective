%Parameters

use_sync = 1;
rec_duration = 3600;                                                                            %approx rec duration (s)
ntw_time = 15.65e-12;                                                                           %network tic interval (s)
n_tags = 5;
TTL_time_diff = [21; 13; 8; 5; 4];                                                              %TTL delays in s
CDPmtdata.Fs = 100;    

data = load('extracted_210222_cdp_1.mat');


for i = 1:1
    bat_pos = data.tag_data{i}(1:350000,3:4);
    disp(size(bat_pos));
    
    bat_acc = data.tag_ac_data{i}(1:350000, 3:4);
    bat_acc(:,1) = bat_acc(:,1) - mean(bat_acc(:,1));
    bat_acc(:,2) = bat_acc(:,2) - mean(bat_acc(:,2));
    sigma_a_x = sqrt(var(bat_acc(:,1)));
    sigma_a_y = sqrt(var(bat_acc(:,2)));
    disp(size(bat_acc));
    
    z = cat(2, bat_pos(:,1), bat_acc(:,1)*0.1, bat_pos(:,2), bat_acc(:,2)*0.1);
    disp(size(z));
    [tag_filt, P] = acc_kalman_filter(z, 1/100, 5, 0.5);
    disp(P);
    
end
disp("filter done!")

samples = linspace(1,length(bat_pos),length(bat_pos));
plot(samples, tag_filt(1,:), samples, data.tag_data_filt{1}(1:350000, 3));
shg;
%plot(tag_filt(1,:))
%plot(data.tag_data_filt{1}(:, 3))

