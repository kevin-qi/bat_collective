function [pos] = extract_position(data)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
N = size(data.tag_data_filt, 2);

for b_num=1:N
    pos{b_num} = data.tag_data_filt{b_num}(10000:355000,3:5);
end

end

