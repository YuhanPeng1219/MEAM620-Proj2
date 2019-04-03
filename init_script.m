% Add additional inputs after sensor if you want to
% Example:
% your_input = 1;
% estimate_vel_handle = @(sensor) estimate_vel(sensor, your_input);
%
% We will only call estimate_vel_handle in the test function.
% Note that thise will only create a function handle, but not run the function

% clear all;
% load('C:\Users\rainc\OneDrive\Desktop\meam 620\proj2phase2\data\studentdata1.mat')
% result1 = [];
% result2 = [];
estimate_vel_handle = @(sensor) estimate_vel(sensor);

% diff = [];
% for a = 1:size(time,2)-1
%     diff = [diff;time(a+1)-time(a)];
% end
% 
% for i = 1:size(data,2)
%    sensor = data(i);
%    
%    [vel, omg] = estimate_vel_handle(sensor);
%    
%    quad_info = [vel;omg];
%    result1 = [result1, quad_info];
% end