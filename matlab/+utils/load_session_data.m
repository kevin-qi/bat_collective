function [output] = load_session_data(sessions)
%LOAD_SESSION_DATA Load data from mutliple sessions into single cell array
%indexed in order of date (oldest to newest)

session_data = {};
disp(size(sessions,1));
for sess_index = 1:size(sessions,1)
    disp(sess_index);
    sess_name = sessions(sess_index, :);
    data_path = sprintf('../data/%s/rtls/Analysis_%s/Analysis_%s.mat', sess_name, sess_name, sess_name)
    disp(data_path);
    data = load(data_path);

    session_data{sess_index} = data;
end

output = session_data;
end

