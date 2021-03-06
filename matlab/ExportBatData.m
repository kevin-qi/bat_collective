%% Export bat data as .csv for python compatibility

%data = load('data/210305/rtls/Analysis_210305/Analysis_210305.mat');

N = 5;

bat_pos_data = zeros([30000, N*3]);
for i = 1:5
    bat_pos_data(:, 3*i-2:3*i) = data.tag_data_filt{i}(20001:50000,3:5);
end

%bflying_data(:, :) = test_data.bflying(1:350000,:);

%writematrix(bat_pos_data, '210305_bat_pos_mat.csv');
%writematrix(bflying_data, '210305_bflying_mat.csv');

writematrix(bat_pos_data, 'test_pos_data.csv');