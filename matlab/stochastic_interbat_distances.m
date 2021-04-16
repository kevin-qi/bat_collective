import utils.*;

figure;
config1 = [[-1 -1];[1 -1];[0 1];[1 1]];
config2 = [[-1 -1];[1 -1];[0 -1]];
pos{1} = datasample(config1, 10000, 1).';
pos{2} = datasample(config2, 10000, 1).';
dist = sqrt(sum((pos{2} - pos{1}).^2,1));
histogram(dist);

config2 = [[-1 -1];[1 -1]; [0 1]];