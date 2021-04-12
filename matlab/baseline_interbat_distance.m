
import voronoi_utils.*;
import utils.*;

%Parameters
x2=2.8; x1=-2.8; y2=2.8;  y1=-2.8;  z1=0; z2=2.30;             %Flight volume coordinates
edges_d = {x1:(x2-x1)/10:x2 y1:(y2-y1)/10:y2};                 %Edges for density histogram
bowl = [-0.29, 0.05, 0.45];                                      %x,y,z bowl 1
Fs = 100;                                                      %resampling frequency (Hz) for common time
n_tags = 5;
bat_nms = ['Dai'; 'Den'; 'Dia'; 'Dor'; 'Dum'];
bat_pairs = nchoosek(1:n_tags,2);   bat_pair_nms = [bat_nms(bat_pairs(:,1),:), '-'.*ones(length(bat_pairs),1), bat_nms(bat_pairs(:,2),:)];
bat_clr = lines(n_tags);
v_th = 0.5;                                                    %Velocity threshold (m/s)
TTL_time_diff = [21; 13; 8; 5; 4];
N = 5;
max_distance = norm([x2-x1 y2-y1]);

N_points = 10000;
bounds = [[linspace(x1,x2,N_points);linspace(y2,y2,N_points)]...
          [linspace(x2,x2,N_points);linspace(y2,y1,N_points)]...
          [linspace(x2,x1,N_points);linspace(y1,y1,N_points)]...
          [linspace(x1,x1,N_points);linspace(y1,y2,N_points)]];
      
bounds = unique(bounds.','rows').';

%bounds = [[linspace(x1,x2,N_points);linspace(y2,y2,N_points)]];
N_samples = 350000;
rand_pos = [];
for b_num=1:6
    rand_pos{b_num} = datasample(bounds, N_samples, 2)';
end

pairwise_dist = pairwise_distance(rand_pos);
disp(size(pairwise_dist));
figure;
axes = [];
index = 1;
bins = linspace(0,max_distance,20);

i=1;
j=1;
ax = subplot(1, 1, index);
axes = [axes ax];
dist = min(pairwise_dist{j}, [], 2);
%dist = pairwise_dist{j}(:,2);
disp(mean(dist));
H = histogram(dist, bins, 'Normalization','probability');

%xline(max_distance, 'LineWidth', 1, 'Color', 'r');
if i == 1
    title(bat_nms(j,:), bat_nms(k,:));
end
set(gca,'xtick',[1,2,3,4,5,6,7,8]);
if index > 1
    set(gca,'ytick',[])
end
if index == 1
    set(gca, 'ytick',[0.1,0.2,0.3,0.4,0.5]);
end
index = index+1;

axis tight;

linkaxes(axes, 'xy');
ylim([0,0.5]);
sgtitle('Pairwise Distance Distribution');
 