function [vel, omg] = estimate_vel(sensor, varargin)
%ESTIMATE_VEL 6DOF velocity estimator
%   sensor - struct stored in provided dataset, fields include:
%          - is_ready: logical, indicates whether sensor data is valid; if the
%                      sensor data is invalid, return empty arrays for vel, omg
%          - t: timestamp
%          - rpy, omg, acc: imu readings (you should not use these in this phase)
%          - img: uint8, 240x376 grayscale image
%          - id: 1xn ids of detected tags; if no tags are present return empty
%                arrays for vel, omg
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
%              estimate_vel_handle = ...
%                  @(sensor) estimate_vel(sensor, your personal input arguments);
%   vel - 3x1 velocity of the quadrotor in world frame
%   omg - 3x1 angular velocity of the quadrotor

vel = zeros(3,0);
omg = zeros(3,0);

persistent First_time img_old t_old V_old delta_old

%% ***************** Determine whether sensor is ready ********************
if sensor.is_ready == 0
    warning("Sensor isn't ready now!!!!!!!!!")
    return
end

%% ********************* Determine whether tag is seen ********************
if isempty(sensor.id) == 1
    t_old = sensor.t;
    img_old = sensor.img;
    warning("No tag is seen!!!!!!!!!!!!")
    return
end

%% ********************** Initialize parameter ****************************
% intrinsic matrix
K = [311.0520, 0, 201.8724;
     0, 311.3885, 113.6210;
     0, 0, 1];     

sample_size = 80;
%% ************************ If the first time  ****************************
if isempty(First_time)
   img_old = sensor.img;
   t_old = sensor.t;
   First_time = 1;
   V_old =[0;0;0;0;0;0];
   delta_old = 0;
   return
end

%% ****************** Detect corners from img *****************************
img_new = sensor.img;
corners = detectFASTFeatures(img_old);%, 'MinQuality', 0.1);       % Detect features
imgPts1 = corners.selectStrongest(sample_size).Location;

%% ************************ Track corner **********************************
% Create tracker
track = vision.PointTracker; %('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);

% initialize tracker with detected corner and initial frame
initialize(track, imgPts1, img_old);

% track conrners
[imgPts2, validIndex] = track.step(img_new);
matchPts1 = imgPts1(validIndex, :);
matchPts2 = imgPts2(validIndex, :);

% Used for debug
% figure;
% showMatchedFeatures(img_old, img_new, matchPts1, matchPts2);

% update img_old
img_old = img_new;

%% ********************* compute optical flow *****************************
% convert to projective space
matchPts1 = [matchPts1,ones(size(matchPts1, 1),1)]';
matchPts2 = [matchPts2,ones(size(matchPts2, 1),1)]';
% imgPts1 = [imgPts1,ones(size(imgPts1, 1),1)]';
% imgPts2 = [imgPts2,ones(size(imgPts2, 1),1)]';

% compute calibrated coordinate
Pts1Cal = K \ matchPts1;            
Pts2Cal = K \ matchPts2;                % 3 by n
% Pts1Cal = K \ imgPts1;            
% Pts2Cal = K \ imgPts2;  

a=0.8;

% time stamp
delta_t = a*(sensor.t - t_old) + (1-a)*delta_old;
% delta_t = sensor.t - t_old;
if delta_t < 0.02
    delta_t = 0.02;
end
t_old = sensor.t;
delta_old = delta_t;

% optical flow
Opt_flow = (Pts2Cal-Pts1Cal)/delta_t;   % 3 by n
% test = abs(Opt_flow) > 3.5;
% [~,b] = find(test);
% Opt_flow(:,b) = []; Pts2Cal(:,b) = [];

%% ********************** Compute depth ***********************************
% get T,R
[T, R, ~, ~] = estimate_pose_new(sensor);

% compute homography
Homography = [R(:,1:2),T];
Homography = Homography / Homography(end, end);

% Real world coordinate
WCoord = (Homography) \ Pts2Cal;  
WCoord = WCoord ./ WCoord(3, :);    % 3 by n
WorldCoord = [WCoord(1:2,:); zeros(1, size(WCoord, 2));WCoord(3,:)]; % 4 by n

% compute the coordinate in camera
World2Cam_homo = [R,T;0,0,0,1];
CameraCoord = World2Cam_homo * WorldCoord;

% Depth
Z = CameraCoord(3, :);

%% *********************** Construct Matrix F *****************************
fmat = [];
for i = 1:size(Pts2Cal, 2)
    fmat = [fmat; compute_f(Z(i), Pts2Cal(1:2, i))];
end

%% ***************** Use Ransac to reject outliers ************************
% reshape the flow 
xydot = Opt_flow(1:2,:);
xydot = xydot(:);
% V = fmat\xydot;

[flow, F_new] =  estimateFmat_Ransac(Pts2Cal, 3, xydot, fmat);

% figure out V
V = F_new\flow;

% beta = 0.632; % template
beta = 0.6; gama = 0.96;
V(1:3) = beta*V(1:3) + (1-beta)*V_old(1:3);
V(4:5) = gama*V(4:5) + (1-gama)*V_old(4:5);

V_old = V; 

vel = V(1:3);
omg = V(4:6);

% % figure out V
% V = fmat\xydot;
% vel = V(1:3);
% omg = V(4:6);

vel = R'*vel;
omg = R'*omg;

end

%% ******************* function to compute f1 f2 **************************
function f = compute_f(Z, pt)
    x = pt(1); y =pt(2);
    A = (1/Z)*([-1,0,x;
             0,-1,y]);
    B = [x*y, -(1+x^2), y; 
         1+y^2, -x*y, -x];
    f = [A,B];
end

%% ******************** Using Ransac to Reject outliers *******************
% Do: use Ransac to reject outliers and construct new F and flow. The
%     outliers are rejected according to the norm of error of flow
%     Input:
%           xy: points to construct f matrix 2 by n
%           sample_size: the minimum sample_size for ransac
%           flow_old: flow waiting to update 2 by n
%           fmat_old: F matrix 2*n by 6
%     Output:
%           flow: the updated flow, less than flow_old
%           F_new: the new matrix F constructed by inliers N by 6 
%
function [flow, F_new] =  estimateFmat_Ransac(xy, sample_size, flow_xy, fmat_old)
    % initialize parameter
    maxIterat = 200;
    sz = sample_size;
    eps = 0.11;
    MostInliers = 0;
    
%     maxIterat = 300;
%     sz = sample_size;
%     eps = 0.05;
%     MostInliers = 0;
    
    % main  loop to randomly sample
    for i = 1:maxIterat
        % randomly sort
        indices = randperm(size(xy,2));
        sampleInd = indices(1:sz);
        testInd =  indices(sz+1:length(indices));  
        
        % compute V using sampleInd
        sapInd = [2*sampleInd-1;2*sampleInd];
        sapInd = sapInd(:);
        sample_flow = flow_xy(sapInd,:);
        sample_f = fmat_old(sapInd,:);
        sample_V = sample_f\sample_flow;
        
        % compute distance using sample_V 
        % The distance here is norm of [flow - test_f*sample_V]
        reqInd = [2*testInd-1;2*testInd];
        reqInd = reqInd(:);
        test_f = fmat_old(reqInd,:);
        test_flow = flow_xy(reqInd,:);
        test_result = test_f*sample_V;
        test_flow = reshape(test_flow, [2,length(testInd)]);
        test_result = reshape(test_result, [2,length(testInd)]);
        error = vecnorm(test_result-test_flow);
        
        % find inliers
        lower_mat = error < eps;
        Inliers = [sampleInd, testInd(lower_mat)];    
        needInd = [2*Inliers-1;2*Inliers];
        needInd = needInd(:);

        NInliers = length(Inliers);

        if (NInliers > MostInliers)
            MostInliers = NInliers;
%             bestInliers = Inliers;
            flow = flow_xy(needInd,:);
            F_new = fmat_old(needInd,:);
           if MostInliers > 65
               break
           end
        end        
    end
end
