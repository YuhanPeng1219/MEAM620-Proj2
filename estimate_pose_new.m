function [ T, R, pos, q] = estimate_pose_new(sensor, varargin)
%ESTIMATE_POSE 6DOF pose estimator based on apriltags
%   sensor - struct stored in provided dataset, fields include:
%          - is_ready: logical, indicates whether sensor data is valid
%          - rpy, omg, acc: imu readings (you should not use these in this phase)
%          - img: uint8, 240x376 grayscale image
%          - id: 1xn ids of detected tags, if no tags are present return empty
%                arrays for pos, q
%          - p0, p1, p2, p3, p4: 2xn pixel position of center and
%                                four corners of detected tags
%            Y
%            ^ P3 == P2
%            | || P0 ||
%            | P4 == P1
%            o---------> X
%   varargin - any variables you wish to pass into the function, could be
%              a data structure to represent the map or camera parameters,
%              your decision. But for the purpose of testing, since we don't
%              know what inputs you will use, you have to specify them in
%              init_script by doing
%              estimate_pose_handle = ...
%                  @(sensor) estimate_pose(sensor, your personal input arguments);
%   pos - 3x1 position of the quadrotor in world frame
%   q   - 4x1 quaternion of the quadrotor [w, x, y, z] where q = w + x*i + y*j + z*k

pos = zeros(3,0);
q = zeros(4,0);

%% ***************** Determine whether sensor is ready ********************
if sensor.is_ready == 0
    warning("Sensor isn't ready now!!!!!!!!!")
    return
end

%% ********************* Determine whether tag is seen ********************
if isempty(sensor.id) == 1
    warning("No tag is seen!!!!!!!!!!!!")
    return
end

%% ************************** Initiate parameter **************************
% intrinsic matrix
K = [311.0520, 0, 201.8724;
     0, 311.3885, 113.6210;
     0, 0, 1];     

% tag distance in real world
dist_l = 0.178;
dist_s = 0.152;


%% ********************** Compute required position ***********************

[row, col] = id2sub(sensor.id);
tagWorldPos = compute_worldPos(row, col, sensor.id, dist_l, dist_s);
tagImgPixel = assign_ImgPixel(sensor);

%% ********************* Construct A matrix *******************************
A = zeros(10*size(sensor.id, 2), 9);
for k = 1:size(tagWorldPos, 3)
    for j = 1:size(tagWorldPos, 2)
        row_num = 2*(k-1)*size(tagWorldPos, 2)+2*j-1;
        A(row_num:row_num+1,:) = Aequ_compute(tagWorldPos(:,j,k), tagImgPixel(:,j,k));
    end
end

%% ******************** Estimate Homography *******************************
[~, ~, V] = svd(A);
H = reshape(V(:,end), [3,3]);
H = H';
H = H / H(end, end);
 
%% ************************ Estimate R,T **********************************
Cal_H = K\H;
R_est = [Cal_H(:, 1),Cal_H(:, 2), cross(Cal_H(:,1),Cal_H(:,2))];
[u,~,v] = svd(R_est);
R = u*[1,0,0;0,1,0;0,0,det(u*v')]*v';
T = Cal_H(:,3)/(0.5*(norm(Cal_H(:,1))+norm(Cal_H(:,2))));

%% ********************** Calculate pos and ori ***************************
%% Calculate pos and q
Homo_trans = [R,T;0,0,0,1];
PosInCam = [-0.04; 0.0; -0.03; 1];
pos = (Homo_trans) \ PosInCam;
pos(4)=[];
RotCamtoQuad = rotx(180)*rotz(44.5);
RotQuadinWorld = (RotCamtoQuad \ R)';
q = (rot2quater(RotQuadinWorld))';
end

%% ********************* Compute sub of sensor id *************************

function [row, col] = id2sub(id)
   col = zeros(1, size(id, 2));
   row = zeros(1, size(id, 2));
   col = ceil((id+1)/12);      % compute col
   row = id+1 - (col-1)*12;    % compute row
end

%% *********************** Compute pos of tag *****************************
function tagWorldPos = compute_worldPos(row, col, id, dl, ds)
    tagWorldPos = zeros(2, 5, size(id, 2));
    % compute P0 first
    for i = 1 : size(id, 2)
        if col(1, i) < 4
            pos0_x = ds/2 + 2*ds*(row(1, i)-1);
            pos0_y = ds/2 + 2*ds*(col(1, i)-1);
            tagWorldPos(:,1,i) = [pos0_x; pos0_y];         
        elseif col(1, i) <7
            pos0_x = ds/2 + 2*(row(1, i) - 1)*ds;
            pos0_y = ds/2 + (2*(col(1, i) - 1)-1)*ds + dl;
            tagWorldPos(:,1,i) = [pos0_x;pos0_y];
        else  
            pos0_x = ds/2 + 2*(row(1,i) - 1)*ds;
            pos0_y = ds/2 + (2*(col(1,i) - 1)-2)*ds + 2*dl;
            tagWorldPos(:,1,i) = [pos0_x;pos0_y];
        end
        % compute P1 P2 P3 P4
            tagWorldPos(:,2,i) = tagWorldPos(:,1,i) + ds/2*[1; -1];
            tagWorldPos(:,3,i) = tagWorldPos(:,1,i) + ds/2*[1; 1];
            tagWorldPos(:,4,i) = tagWorldPos(:,1,i) + ds/2*[-1; 1];
            tagWorldPos(:,5,i) = tagWorldPos(:,1,i) + ds/2*[-1; -1];  
    end
end

%% ********************* Assign the img pixel *****************************
function tagImgPixel = assign_ImgPixel(sensor)
   tagImgPixel = zeros(2, 5, size(sensor.id, 2));
   tagImgPixel(:,1,:) = sensor.p0;
   tagImgPixel(:,2,:) = sensor.p1;
   tagImgPixel(:,3,:) = sensor.p2;
   tagImgPixel(:,4,:) = sensor.p3;
   tagImgPixel(:,5,:) = sensor.p4;
end

%% ********************* Function to compute A equation *******************
function A_equation = Aequ_compute(X, Y)
    x = X(1);
    y = X(2);
    x_o = Y(1);
    y_o = Y(2);
    A_equation = [x, y, 1, 0, 0, 0 ,-x*x_o, -y*x_o, -x_o;...
                  0, 0, 0, x, y, 1, -x*y_o, -y*y_o, -y_o];
end

%% ****************** Function to compute quarternion *********************
function q = rot2quater(R)
q = zeros(4,0);
tau = trace(R);
cosphi = 0.5*(tau-1);
phi = acos(cosphi);
u_hat = (R-R')/(2*sin(phi));
u = solve_hat(u_hat);
q = [cos(phi/2), sin(phi/2)*u'];
end

function t=solve_hat(t_hat)
    t=[t_hat(3,2); t_hat(1,3); t_hat(2,1)];
end
