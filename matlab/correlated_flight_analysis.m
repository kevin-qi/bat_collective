%% Script for the analysis of collective flight behavior

%data = load('data/Analysis_210222/Analysis_210222.mat');
bflying = data.bflying;
N = min(size(bflying));
fs = 100;
total_flying = sum(bflying, 2);
%plot(total_flying);
%shg;

time = linspace(1, length(total_flying)/100, length(total_flying));
T = length(time);
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

windows = reshape(x_pos

