clear all;
close all;
clc;

author = "Gabriel Gekas";

heading = 180; % degrees heading in geodetic frame (north-0, +clockwise)
tilt = 30; % degrees upward tilt of hydrophone array

% water properties taken from https://data.caloos.org/
T = 18.39; % Temperature [Celsius]
S = 30.2; % Salinity [PSS]
c = 1500; % m/s

array_axis_size = 2; % m

% define hydrophone placement in array coordinates [meters]
A = [0.45 0.54 0.11]';
B = [0.44 0.03 0.61]';
C = [0.00 0.44 0.39]';
D = [0.43 0.44 0.60]';
E = [0.60 0.02 0.02]';

phone_coords = [A, B, C, D, E]; % define matrix of coordinates

% dir = [cosd(90 - heading), sind(90 - heading), 0];
dir = [1, -1, 0];
dir = dir / norm(dir);
xang = 180;
yang = 0;
zang = 135 - heading;
R = rodrigues(dir, tilt) * eul2rotm([deg2rad(zang), deg2rad(yang), deg2rad(xang)], "ZYX");
p = zeros(3, 1); % get position of array
phones = (pinv(R) * phone_coords)';

x = phones(:, 1);
y = phones(:, 2);
z = phones(:, 3);

% figure; 
% plot3(phone_coords(1, :), phone_coords(2, :), phone_coords(3, :), 'o');
% xlabel('X')
% ylabel('Y')
% zlabel('Z')

%% calculate example TDOA:
% comment this section out if TDOA is supplied elsewhere

source = [1, -6, -4]'; % 10 meters south
r = phones' - source;
r = sqrt(r(1, :).^2 + r(2, :).^2 + r(3, :).^2); % calculates range to each receiver
t = r/c; % get travel time to each receiver
t = t - t(1); % get reduced time relative to hydrophone A
tA = t(1); tB = t(2); tC = t(3); tD = t(4); tE = t(5);


%% calculate position from from TDOA:

tdoaest = t(2:end); % tdoas referenced to t(1) (K)x(L-1)
tdoavar = 0.00001 * ones(1, size(tdoaest, 2));
anchorpos = phones'; % (Q)x(L) where Q = 3 dimensions

[tgtposest, tgtposcov] = tdoaposest(tdoaest, tdoavar, anchorpos, PropagationSpeed = c);
det = tgtposest';
detvar = 5e-7 * tgtposcov;

% get bearing and elevation angle
[azimuth, elevation, r] = cart2sph(det(1), det(2), det(3));
bearing = 90 - rad2deg(azimuth);
elevation = rad2deg(elevation);

% p1 = inv(R)*B;
% p2 = inv(R)*C;
% dt = tB - tC;
% [ang, det] = TDOA2azimuth(p1, p2, dt, c);
% det = (R*det')' % detection in array frame


%% plot array and reception

x_all = [x; det(1)];
y_all = [y; det(2)];
z_all = [z; det(3)];

width = 1000; height = 400;
Pix_SS = get(0,'screensize');
fig = figure('Position', [(Pix_SS(3)-width)/2 (Pix_SS(4)-height)/2 width height]); a = tiledlayout(3, 7);
axis_offset = 2; % m
lat_offset = 2; % m
set(gcf, 'Color', [1, 1, 1]);
title(a, 'Acoustic Array Pose in World Frame (ENU frame)', 'Interpreter', 'Latex');

nexttile([3, 3]);
tp1 = theaterPlot('Parent', fig.CurrentAxes,...
    'XLim',[min(x_all)-lat_offset max(x_all)+lat_offset],...
    'YLim',[min(y_all)-lat_offset max(y_all)+lat_offset],...
    'ZLim',[min(z_all)-axis_offset max(z_all)+axis_offset]);
op1 = orientationPlotter(tp1,'DisplayName','array','LocalAxesLength', array_axis_size);
lp1 = detectionPlotter(tp1,'DisplayName','hydrophones', 'VelocityScaling', 10);
det1 = detectionPlotter(tp1,'DisplayName','detection', 'VelocityScaling', 10, 'MarkerFaceColor', 'b');
set(fig.CurrentAxes, 'Xdir', 'normal');
set(fig.CurrentAxes, 'Ydir', 'normal');
set(fig.CurrentAxes, 'Zdir', 'normal');
view(0, 90) % views in line with z axis

nexttile([3, 2]);
tp2 = theaterPlot('Parent', fig.CurrentAxes,...
    'XLim',[min(x_all)-lat_offset max(x_all)+lat_offset],...
    'YLim',[min(y_all)-lat_offset max(y_all)+lat_offset],...
    'ZLim',[min(z_all)-axis_offset max(z_all)+axis_offset]);
op2 = orientationPlotter(tp2,'DisplayName','array','LocalAxesLength', array_axis_size);
lp2 = detectionPlotter(tp2,'DisplayName', 'hydrophones', 'VelocityScaling', 10);
det2 = detectionPlotter(tp2,'DisplayName','detection', 'VelocityScaling', 10, 'MarkerFaceColor', 'b');
set(fig.CurrentAxes, 'Xdir', 'normal');
set(fig.CurrentAxes, 'Ydir', 'normal');
set(fig.CurrentAxes, 'Zdir', 'normal');
view(0, 0) % views in line with x axis

nexttile([3, 2]);
tp3 = theaterPlot('Parent', fig.CurrentAxes,...
    'XLim',[min(x_all)-lat_offset max(x_all)+lat_offset],...
    'YLim',[min(y_all)-lat_offset max(y_all)+lat_offset],...
    'ZLim',[min(z_all)-axis_offset max(z_all)+axis_offset]);
op3 = orientationPlotter(tp3,'DisplayName','array','LocalAxesLength', array_axis_size);
lp3 = detectionPlotter(tp3,'DisplayName','hydrophones', 'VelocityScaling', 10);
det3 = detectionPlotter(tp3,'DisplayName','detection', 'VelocityScaling', 10, 'MarkerFaceColor', 'b');
set(fig.CurrentAxes, 'Xdir', 'normal');
set(fig.CurrentAxes, 'Ydir', 'normal');
set(fig.CurrentAxes, 'Zdir', 'normal');
view(90, 0) % views in line with y axis


% orientation and position
plotOrientation(op1, R, p');
plotDetection(lp1, phones);
plotDetection(det1, det, detvar);

% plotDetection(det1, p1');
% plotDetection(det1, p2');

% orientation and position
plotOrientation(op2, R, p');
plotDetection(lp2, phones);
plotDetection(det2, det, detvar);

% orientation and position
plotOrientation(op3, R, p');
plotDetection(lp3, phones);
plotDetection(det3, det, detvar);



%% functions

function inv_T = invT(T)
    % calculates the inverse of a homogenous tranformation matrix
    top_left = T(1:3, 1:3)';
    top_right = -top_left*T(1:3, 4);
    bottom = [0, 0, 0, 1];
    inv_T = vertcat(horzcat(top_left, top_right), bottom);
end

function [R, p] = T2RP(T)
    R = T(1:3, 1:3);
    p = T(1:3, 4);
end

function rod_mat = rodrigues(x, theta)
    % given a rotation axis and angle, calculates rotation matrix
    if length(x) ~= 3
        sprintf("Input must be a 3x1 or 1x3 vector")
        return
    end
    rod_mat = eye(3) + sind(theta)*skew(x) + (1 - cosd(theta))*skew(x)^2;
end

function skew_mat = skew(x)
    % given a 3x1 position matrix, the skew matrix is calculated
    if length(x) ~= 3
        sprintf("Input must be a 3x1 or 1x3 vector")
        return
    end
    skew_mat = [0,    -x(3),  x(2); ...
                x(3),     0, -x(1); ...
               -x(2),  x(1),     0];
end

function [ang, xyz] = TDOA2azimuth(p1, p2, dt, c)
% computes the angle made between the ray path and the axis between the
% hydrophones
% :param p1: position of hydrophone 1
% :param p2: position of hydrophone 2
% :param dt: time difference between arrivals (t2 - t1)
% :param c: sound speed (m/s)

% sqrt((xt-x2)^2 + (yt-y2)^2 + (zt-z2)^2) -  sqrt((xt-x1)^2 + (yt-y1)^2 +
% (zt-z1)^2) = c*(t2-t1)
mid = (p1+p2)./2; % calculate midpoint of p1, p2

syms x y
z = 0;

eqn = sqrt((x - p2(1))^2 + (y - p2(2))^2 + (z - p2(3))^2)...
    - sqrt((x - p1(1))^2 + (y - p1(2))^2 + (z - p1(3))^2)...
    == c*dt;
guess = [0; -10]; % hard coded to search south
xy = vpasolve(eqn, [x, y], guess);
y = double(xy.y);
x = double(xy.x);
xyz = [x, y, z];

Y = y-mid(2);
X = x-mid(1);

ax = atan2d(p2(2) - p1(2), p2(1) - p1(1)); % axis in xy plane that the hydrophones lie on
% ad = 180 - ax;
ang = 90 - atan2d(Y, X); % heading of detection
% ang = mod(ang+180, 360); 

end