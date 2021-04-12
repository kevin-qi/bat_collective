function [rand_pos] = rand_bat_bound_pos(N_samples, bounds)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x1 = bounds.x1;
x2 = bounds.x2;
y1 = bounds.y1;
y2 = bounds.y2;

z_pos = 2;

N_points = 1000;
bounds = [[linspace(x1,x2,N_points);linspace(y2,y2,N_points)]...
          [linspace(x2,x2,N_points);linspace(y2,y1,N_points)]...
          [linspace(x2,x1,N_points);linspace(y1,y1,N_points)]...
          [linspace(x1,x1,N_points);linspace(y1,y2,N_points)]];
      
bounds = unique(bounds.','rows').';

rand_pos = datasample(bounds, N_samples, 2)';
rand_pos = [rand_pos linspace(z_pos,z_pos,N_samples)'];
end

