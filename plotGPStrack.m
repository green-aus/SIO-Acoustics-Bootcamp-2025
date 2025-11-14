clear;
clc;
close all;


%% read gps data
gps1file = 'Data/SP2423_gnss_gp170-2024-11-02.txt';
gps2file = 'Data/SP2423_gnss_gp170-2024-11-03.txt';

% gps1file = 'C:\Users\gabe\Desktop\School\PhD\FA24\SIOB 280\cruise_report\SP2424\openrvdas\data\SP2424_gnss_gp170-2024-11-03.txt';
% gps2file = 'C:\Users\gabe\Desktop\School\PhD\FA24\SIOB 280\cruise_report\SP2424\openrvdas\data\SP2424_gnss_gp170-2024-11-04.txt';

titlestring = 'Cruise Track with Bathymetry 11/02/2024';


%read first file
[date,position] = textread(gps1file,'%s %s');
%find indices with $GPGGA
rowIndices = find(contains(position, 'GPGGA'));
pos1 = position(rowIndices);
date1 = date(rowIndices);

%read second file
[date,position] = textread(gps2file,'%s %s');
%find indices with $GPGGA
rowIndices = find(contains(position, 'GPGGA'));
pos2 = position(rowIndices);
date2 = date(rowIndices);

allpos = [pos1;pos2];
alldates = [date1;date2];

for a = 1:length(allpos)
    parts = strsplit(allpos{a}, ',');
    % Assign the individual parts to variables
    [ID, time{a,1}, latitude{a,1}, lat, longitude{a,1}, lon, v1{a,1}, v2{a,1}, v3{a,1}, ...
        v4{a,1}, M, v5{a,1}, M2, other{a,1}] = deal(parts{:});
end

latitude = cellfun(@str2double, latitude);
longitude = cellfun(@str2double, longitude);

% format in the data is 32 degrees 38.3352 minutes N and 117 degrees, 31.2612 minutes W
% it never goes below or above 32, so save to hard code this solution
latMinutes = latitude - 3200;
lonMinutes = longitude - 11700;
decLatitude = 32 + latMinutes/60; %final cruise track latitudes
decLongitude = 117 + lonMinutes/60;
decLongitude = decLongitude * (-1); % 117 W -> -117 final cruise track longitudes

% convert ISO date into Matlab serial date
alldates = dbISO8601toSerialDate(alldates);

%% plot map
% Load bathymetry data
ncfilename = [pwd, '/Gebco_bathymetry-20241119T171750Z-001/Gebco_bathymetry/gebco_2024_n35.3512_s31.7079_w-121.7881_e-116.8022.nc'];
lon_bath = ncread(ncfilename, 'lon'); % Longitude array
lat_bath = ncread(ncfilename, 'lat'); % Latitude array
depth_bath = ncread(ncfilename, 'elevation'); % Depth array
depth_bath = double(depth_bath.');
depth_bath(depth_bath == 32767) = NaN;  % Replace the missing data with NaN

% Convert longitude and latitude to a grid for plotting
[lonGrid, latGrid] = meshgrid(lon_bath, lat_bath);

% Define map limits
latlim = [32.6, 33.0]; % Adjust for your region of interest
lonlim = [-117.5 -117.1];

% Plot bathymetry map
width = 800; height = 500;
Pix_SS = get(0,'screensize');
fig = figure('Position', [(Pix_SS(3)-width)/2 (Pix_SS(4)-height)/2 width height]); a = tiledlayout(3, 7);
set(gcf, 'Color', [1, 1, 1]);
axesm('Mercator', 'Frame', 'off', 'Grid', 'on', 'ParallelLabel', 'off');  % Set Mercator projection
hold on
worldmap(latlim, lonlim)
geoshow(latGrid, lonGrid, depth_bath,'DisplayType','surface');

% define colormap
zlimits = [-1500 500];
demcmap(zlimits)
colorbar

% Overlay cruise track
scatterm(decLatitude, decLongitude, 5, alldates, 'filled','r')

% Label plot
title(titlestring)
xlabel('Longitude')
ylabel('Latitude')
hold off

DIRIN = [pwd, '/figures/'];
filename = 'cruise_track.png';
filepath = [DIRIN, filename];
saveas(gcf, filepath)