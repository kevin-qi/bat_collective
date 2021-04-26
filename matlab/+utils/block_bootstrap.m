function [samples] = block_bootstrap(timeseries, block_size, num_samples)
%BLOCK_BOOTSTRAP Summary of this function goes here
%   Detailed explanation goes here
num_blocks = floor(length(timeseries) / block_size);

chunked_ts = reshape(timeseries(1:block_size*num_blocks), [block_size num_blocks]).';
samples = datasample(chunked_ts, num_blocks * num_samples, 1);
%samples = reshape(samples.', [block_size num_blocks num_samples]);
samples = reshape(samples.', [block_size*num_blocks num_samples]).';


end

