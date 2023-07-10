

% Uncomment/comment lines below to choose time resolution

% file_1 = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSBard_1min';
% file_2 = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSBerg_1min';
% file(2).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSPolargirl_1min';
% file(3).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSBillefjord_1min';
% file(4).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Narveneset_1min';
% file(5).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Gasoyane_1min';
% file(6).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Daudmannsodden_1min';
% file(7).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Bohemanneset_1min';


file_1 = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSBard_10min';
file_2 = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSBerg_10min';
file(2).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSPolargirl_10min';
file(3).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSBillefjord_10min';
file(4).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Narveneset_10min';
file(5).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Gasoyane_10min';
file(6).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Daudmannsodden_10min';
file(7).name = 'https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Bohemanneset_10min';

for id = 1:7

    if id == 1
            a = ncread(file_1,'time')./(24*60*60);          b = ncread(file_2,'time')./(24*60*60);          IWIN(id).time = [a(:)' b(:)']' + datenum('1970-01-01 00:00:00');
            a = ncread(file_1,'temperature');               b = ncread(file_2,'temperature');               IWIN(id).temperature = [a(:)' b(:)']';
            a = ncread(file_1,'relative_humidity');         b = ncread(file_2,'relative_humidity');         IWIN(id).relative_humidity = [a(:)' b(:)']';
            a = ncread(file_1,'air_pressure');              b = ncread(file_2,'air_pressure');              IWIN(id).air_pressure = [a(:)' b(:)']';
            a = ncread(file_1,'longitude');                 b = ncread(file_2,'longitude');                 IWIN(id).lon = [a(:)' b(:)']';
            a = ncread(file_1,'latitude');                  b = ncread(file_2,'latitude');                  IWIN(id).lat = [a(:)' b(:)']';
            a = ncread(file_1,'GPS_speed');                 b = ncread(file_2,'GPS_speed');                 IWIN(id).GPS_speed = [a(:)' b(:)']';
            a = ncread(file_1,'GPS_heading');               b = ncread(file_2,'GPS_heading');               IWIN(id).GPS_heading = [a(:)' b(:)']';
            a = ncread(file_1,'exhaust_plume_influence');   b = ncread(file_2,'exhaust_plume_influence');   IWIN(id).exhaust_plume_influence = [a(:)' b(:)']';
            a = ncread(file_1,'wind_speed_corrected');      b = ncread(file_2,'wind_speed_corrected');      IWIN(id).wind_speed_corrected = [a(:)' b(:)']';
            a = ncread(file_1,'wind_direction_corrected');  b = ncread(file_2,'wind_direction_corrected');  IWIN(id).wind_direction_corrected = [a(:)' b(:)']';
    elseif id > 1 && id < 4
            IWIN(id).time = ncread(file(id).name,'time')./(24*60*60) + datenum('1970-01-01 00:00:00');
            IWIN(id).exhaust_plume_influence = ncread(file(id).name,'exhaust_plume_influence');
            IWIN(id).wind_speed_corrected = ncread(file(id).name,'wind_speed_corrected');
            IWIN(id).wind_direction_corrected = ncread(file(id).name,'wind_direction_corrected');
            IWIN(id).lon = ncread(file(id).name,'longitude');
            IWIN(id).lat = ncread(file(id).name,'latitude');
            IWIN(id).GPS_speed = ncread(file(id).name,'GPS_speed');
            IWIN(id).GPS_heading = ncread(file(id).name,'GPS_heading');
            IWIN(id).temperature = ncread(file(id).name,'temperature');
            IWIN(id).relative_humidity = ncread(file(id).name,'relative_humidity');
            IWIN(id).air_pressure = ncread(file(id).name,'air_pressure');
    else
            IWIN(id).time = ncread(file(id).name,'time')./(24*60*60) + datenum('1970-01-01 00:00:00');
            IWIN(id).temperature = ncread(file(id).name,'temperature');
            IWIN(id).relative_humidity = ncread(file(id).name,'relative_humidity');
            IWIN(id).air_pressure = ncread(file(id).name,'air_pressure');
            IWIN(id).wind_speed = ncread(file(id).name,'wind_speed');
            IWIN(id).wind_direction = ncread(file(id).name,'wind_direction');
    end

end


% reducing pressure to sea level
IWIN(1).hgt = 20;
IWIN(2).hgt = 20;
IWIN(3).hgt = 20;
IWIN(4).hgt = 7;
IWIN(5).hgt = 19;
IWIN(6).hgt = 39;
IWIN(7).hgt = 12;

for id = 1:length(IWIN)
    IWIN(id).mslp = IWIN(id).air_pressure + IWIN(id).hgt./8;
end


IWIN(1).name = 'Bard';
IWIN(2).name = 'Polargirl';
IWIN(3).name = 'Billefjord';
IWIN(4).name = 'Narveneset';
IWIN(5).name = 'Gasoyane';
IWIN(6).name = 'Daudmannsodden';
IWIN(7).name = 'Bohemanneset';




