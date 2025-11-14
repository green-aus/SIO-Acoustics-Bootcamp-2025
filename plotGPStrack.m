clear;
clc;
close all;

% 
% %% read gps data
% gps1file = 'Data/SP2423_gnss_gp170-2024-11-02.txt';
% gps2file = 'Data/SP2423_gnss_gp170-2024-11-03.txt';
% 
% % gps1file = 'C:\Users\gabe\Desktop\School\PhD\FA24\SIOB 280\cruise_report\SP2424\openrvdas\data\SP2424_gnss_gp170-2024-11-03.txt';
% % gps2file = 'C:\Users\gabe\Desktop\School\PhD\FA24\SIOB 280\cruise_report\SP2424\openrvdas\data\SP2424_gnss_gp170-2024-11-04.txt';
% 
% titlestring = 'Cruise Track with Bathymetry 11/02/2024';
% 
% 
% %read first file
% [date,position] = textread(gps1file,'%s %s');
% %find indices with $GPGGA
% rowIndices = find(contains(position, 'GPGGA'));
% pos1 = position(rowIndices);
% date1 = date(rowIndices);
% 
% %read second file
% [date,position] = textread(gps2file,'%s %s');
% %find indices with $GPGGA
% rowIndices = find(contains(position, 'GPGGA'));
% pos2 = position(rowIndices);
% date2 = date(rowIndices);
% 
% allpos = [pos1;pos2];
% alldates = [date1;date2];
% 
% for a = 1:length(allpos)
%     parts = strsplit(allpos{a}, ',');
%     % Assign the individual parts to variables
%     [ID, time{a,1}, latitude{a,1}, lat, longitude{a,1}, lon, v1{a,1}, v2{a,1}, v3{a,1}, ...
%         v4{a,1}, M, v5{a,1}, M2, other{a,1}] = deal(parts{:});
% end
% 
% latitude = cellfun(@str2double, latitude);
% longitude = cellfun(@str2double, longitude);
% 
% % format in the data is 32 degrees 38.3352 minutes N and 117 degrees, 31.2612 minutes W
% % it never goes below or above 32, so save to hard code this solution
% latMinutes = latitude - 3200;
% lonMinutes = longitude - 11700;
% decLatitude = 32 + latMinutes/60; %final cruise track latitudes
% decLongitude = 117 + lonMinutes/60;
% decLongitude = decLongitude * (-1); % 117 W -> -117 final cruise track longitudes
% 
% % convert ISO date into Matlab serial date
% alldates = dbISO8601toSerialDate(alldates);

%% get Tx data
Tx = readtable([pwd '/Data/Tx.csv']);
lat = Tx.GPS_lat_N;
lon = Tx.GPS_long_W;
coords = rmmissing(table(lat, lon));
coords = unique(coords, 'rows');
coords.lon = -coords.lon;

% coord_array = table2array(coords)';
% 
% cutoff_dist = 0.05;
% for j = 1:numel(coords.lat)
%     sum(coord_array(:, ~j).^2) - sum(coord_array(:, j).^2)
%     if any(sum(coord_array(:, ~j).^2) - sum(coord_array(:, j).^2) < cutoff_dist)
% 
%     end
% end
% coords.lat and coords.lon = column vectors
lat = coords.lat(:);
lon = coords.lon(:);

minDist = 100;   % meters (choose your threshold)

keep = true(size(lat));

for i = 1:length(lat)
    if ~keep(i)
        continue
    end

    for j = i+1:length(lat)
        if keep(j)
            d = distance(lat(i), lon(i), lat(j), lon(j), wgs84Ellipsoid);
            if d < minDist
                keep(j) = false;   % remove point too close
            end
        end
    end
end

lat = lat(keep);
lon = lon(keep);

station = [7, 1, 5, 4, 3, 6, 4]';
coords = table(station, lat, lon);

array_coords = [32.8670661, -117.2574346]'; % pier


%% plot map
% Load bathymetry data
% ncfilename = [pwd, '/Gebco_bathymetry-20241119T171750Z-001/Gebco_bathymetry/gebco_2024_n35.3512_s31.7079_w-121.7881_e-116.8022.nc'];
% lon_bath = ncread(ncfilename, 'lon'); % Longitude array
% lat_bath = ncread(ncfilename, 'lat'); % Latitude array
% depth_bath = ncread(ncfilename, 'elevation'); % Depth array
% depth_bath = double(depth_bath.');
% depth_bath(depth_bath == 32767) = NaN;  % Replace the missing data with NaN

bathy = load('Data/sd1p5secf_v4.mat');

% define title
titlestring = 'SIO Pier and Waypoint Geometry';

% Convert longitude and latitude to a grid for plotting
[lonGrid, latGrid] = meshgrid(bathy.lon, bathy.lat);

% Define map limits
% latlim = [32.6, 33.0]; % Adjust for your region of interest
% lonlim = [-117.5 -117.1];
latlim = [32.86, 32.875]; % Adjust for your region of interest
lonlim = [-117.28 -117.25];

% Plot bathymetry map
width = 800; height = 500;
Pix_SS = get(0,'screensize');
fig = figure('Position', [(Pix_SS(3)-width)/2 (Pix_SS(4)-height)/2 width height]); a = tiledlayout(3, 7);
set(gcf, 'Color', [1, 1, 1]);
axesm('Mercator', 'Frame', 'off', 'Grid', 'on', 'ParallelLabel', 'off');  % Set Mercator projection
hold on
worldmap(latlim, lonlim)
% geoshow(latGrid, lonGrid, depth_bath,'DisplayType','surface');
bathymetry = geoshow(latGrid, lonGrid, -bathy.bathy,'DisplayType','surface');

% define colormap
zlimits = [-600 100];
demcmap(zlimits)
c = colorbar;
c.Label.String = 'depth [m]';
c.Label.Interpreter = 'Latex';
c.Label.VerticalAlignment = 'middle';
c.Label.Position = [0.5, 115];
c.Label.Rotation = 0;

% Overlay cruise track
% scatterm(decLatitude, decLongitude, 5, alldates, 'filled','r')

for j = 1:numel(coords.lat)
    ray = plotm([coords.lat(j), array_coords(1)]', [coords.lon(j), array_coords(2)]', 'LineWidth', 1.2, 'Color', 'm');
    coordinate2 = [coords.lat(j), coords.lon(j)]';
    [bearing(j), range(j)] = coordinates2HeadingDistance(array_coords, coordinate2);
    textm(coords.lat(j), coords.lon(j),...
        {['Bearing: ', num2str(round(bearing(j), 2)), ' deg'],...
        ['Range: ', num2str(round(range(j), 2)), ' m']...
        ['WP', num2str(coords.station(j))]}, 'HorizontalAlignment', 'right')
end
bearing = bearing';
range = range';

% plot Tx points
tx = scatterm(coords.lat, coords.lon, 15, 'filled', 'r');
array = scatterm(array_coords(1), array_coords(2), 30, 'filled', 'k');


% Label plot
legend([bathymetry, array, tx, ray], {'Bathymetry', 'Array', 'Transmission', 'Ray Path'});
title(titlestring, 'Interpreter', 'Latex')
xlabel('Longitude')
ylabel('Latitude')
hold off

DIRIN = [pwd, '/figures/'];
filename = 'cruise_track.png';
filepath = [DIRIN, filename];
saveas(gcf, filepath)

tab = table(station, bearing, range)


%% functions

function [heading, x] = coordinates2HeadingDistance(coordinate1, coordinate2)
% This function returns the heading between two coordinates
% 
% :param coordinate1: first waypoint
% :param coordinate2: second waypoint
% Author: Gabriel Gekas

% calculate m/deg arc distance for lat/lon
E = wgs84Ellipsoid;
lat_deg = distance(coordinate1(1)-0.5,coordinate1(2),coordinate1(1)+0.5,coordinate1(2),E); % m/deg arc length
lon_deg = distance(coordinate1(1),coordinate1(2)-0.5,coordinate1(1),coordinate1(2)+0.5,E); % m/deg arc length

lon_skew = lon_deg / lat_deg;

coordinates = remap_coords(coordinate2, -coordinate1, 0);
coordinates(2) = coordinates(2) * lon_skew;
heading = atan2d(coordinates(1), coordinates(2));
heading = remap_heading(heading);

x = distance(coordinate1(1), coordinate1(2), coordinate2(1), coordinate2(2), E);

end


function xy = remap_coords(coords, origin, theta)
    % coords : 2xn vector
    % origin : the amount of translation, 2x1 vector
    % theta : angle in degrees
    R = [cosd(theta), sind(theta), origin(1);...
        -sind(theta), cosd(theta), origin(2);...
                   0,           0,        1];
    pose = vertcat(coords, ones(1, size(coords, 2)));
    rot_coords = R * pose;
    xy = rot_coords(1:2, :);
end

function h_out = remap_heading(h_in)
    % remaps trigonomtric to geodetic
    h_out = mod(90 - h_in, 360);
end