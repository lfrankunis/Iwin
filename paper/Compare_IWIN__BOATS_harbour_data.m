



%% Matching boat locations on Isfjorden

clc

clear Temperature_diff Humidity_diff Windspeed_diff Winddirection_diff




lonlim = [15.6056 15.6127];
latlim = [78.2279 78.2287];



for id = 1:3


    a = zeros(length(IWIN(id).lon),1); b=a; c=a; d=a;


    ida = datenum([2020 06 18 00 00 00]);
    idb = datenum([2023 06 23 00 00 00]);

    a(IWIN(id).time >= ida & IWIN(id).time <= idb)  = 1;


        b(IWIN(id).lon < 15.6127 & IWIN(id).lon > 15.6056 ) = 1;
        c(IWIN(id).lat > 78.2279 & IWIN(id).lat < 78.2287) = 1;
        d(IWIN(id).GPS_speed < 0.1) = 1;


    IWIN(id).temperature(IWIN(id).exhaust_plume_influence == 1) = nan; 


    dd = a + b + c + d;


    ids(id).n = find(dd == 4);
end


k1 = intersect(IWIN(1).time(ids(1).n),IWIN(2).time(ids(2).n));
k2 = intersect(k1,IWIN(3).time(ids(3).n));

for id = 1:3
    [a b c] = intersect(k2,IWIN(id).time);
    IWINs(id).time = IWIN(id).time(c);
    IWINs(id).lon = IWIN(id).lon(c);
    IWINs(id).lat = IWIN(id).lat(c);
    IWINs(id).temperature = IWIN(id).temperature(c);
    IWINs(id).relative_humidity = IWIN(id).relative_humidity(c);
    IWINs(id).air_pressure = IWIN(id).air_pressure(c);
    IWINs(id).mslp = IWIN(id).mslp(c);
    IWINs(id).wind_speed_corrected = IWIN(id).wind_speed_corrected(c);
    IWINs(id).wind_direction_corrected = IWIN(id).wind_direction_corrected(c);
end


for ii = 1:3


    if ii == 1
        idf = [1 2];
    elseif ii == 2
        idf = [1 3];
    elseif ii == 3
        idf = [2 3];
    end

  mp(ii).idf = idf;
  mp(ii).lon_1 = IWINs(idf(1)).lon; 
  mp(ii).lon_2 = IWINs(idf(2)).lon;
  mp(ii).lat_1 = IWINs(idf(1)).lat; 
  mp(ii).lat_2 = IWINs(idf(2)).lat;

  mp(ii).time_1 = IWINs(idf(1)).time;
  mp(ii).time_2 = IWINs(idf(2)).time;

  mp(ii).temperature_1 = IWINs(idf(1)).temperature;
  mp(ii).temperature_2 = IWINs(idf(2)).temperature;
  mp(ii).relative_humidity_1 = IWINs(idf(1)).relative_humidity;
  mp(ii).relative_humidity_2 = IWINs(idf(2)).relative_humidity;
  mp(ii).mslp_1 = IWINs(idf(1)).mslp;
  mp(ii).mslp_2 = IWINs(idf(2)).mslp;
  mp(ii).wind_speed_corrected_1 = IWINs(idf(1)).wind_speed_corrected;
  mp(ii).wind_speed_corrected_2 = IWINs(idf(2)).wind_speed_corrected;
  mp(ii).wind_direction_corrected_1 = IWINs(idf(1)).wind_direction_corrected;
  mp(ii).wind_direction_corrected_2 = IWINs(idf(2)).wind_direction_corrected;

  mp(ii).name_1 = IWIN(idf(1)).name;
  mp(ii).name_2 = IWIN(idf(2)).name;

    temp_diff(ii) = nanmean(IWINs(idf(1)).temperature - IWINs(idf(2)).temperature);
    rh_diff(ii) = nanmean(IWINs(idf(1)).relative_humidity - IWINs(idf(2)).relative_humidity);
    mslp_diff(ii) = nanmean(IWINs(idf(1)).mslp - IWINs(idf(2)).mslp);
    ws_diff(ii) = nanmean(IWINs(idf(1)).wind_speed_corrected - IWINs(idf(2)).wind_speed_corrected);


    temp_mae(ii) = nanmean(abs(IWINs(idf(1)).temperature - IWINs(idf(2)).temperature));
    rh_mae(ii) = nanmean(abs(IWINs(idf(1)).relative_humidity - IWINs(idf(2)).relative_humidity));
    mslp_mae(ii) = nanmean(abs(IWINs(idf(1)).mslp - IWINs(idf(2)).mslp));
    ws_mae(ii) = nanmean(abs(IWINs(idf(1)).wind_speed_corrected - IWINs(idf(2)).wind_speed_corrected));


a11 = IWINs(idf(1)).wind_direction_corrected;
b11 = IWINs(idf(2)).wind_direction_corrected;

        a = rad2deg(angdiff(deg2rad(a11),deg2rad(b11)));
        a = a(isfinite(a));
            a_abs = abs(a);
        wd_diff(ii) = rad2deg(circ_mean(deg2rad(a),[],1));
        wd_mae(ii) = rad2deg(circ_mean(deg2rad(a_abs),[],1));


    Data_points(ii) = length(mp(ii).relative_humidity_1);

end


    temp_diffs = round(temp_diff,2);
    rh_diffs   = round(rh_diff,2);
%     p_diffs    = round(p_diff,2);
    mslp_diffs = round(mslp_diff,2);
    ws_diffs   = round(ws_diff,2);
%     wd_diffs = round(wd_diff,2);

    temp_maes = round(temp_mae,2);
    rh_maes   = round(rh_mae,2);
%     p_maes    = round(p_mae,2);
    mslp_maes = round(mslp_mae,2);
    ws_maes   = round(ws_mae,2);


%%
clc
clear T
Comparisons = ["Bard vs Polargirl","Bard vs Billefjord","Polargirl vs Billefjord"];

  T = table(Comparisons',temp_diffs',rh_diffs',mslp_diffs',ws_diffs',Data_points');
    T.Properties.VariableNames = ["Comparison","Temperature_diff","RelHumidity_diff","Mslp_diff","Windspeed_diff","Data_points"];

    T

%%
clc
    clear TT
Comparisons = ["Bard vs Polargirl","Bard vs Billefjord","Polargirl vs Billefjord"];

  TT = table(Comparisons',temp_maes',rh_maes',mslp_maes',ws_maes',Data_points');
    TT.Properties.VariableNames = ["Comparison","Temperature_mae","RelHumidity_mae","Mslp_mae","Windspeed_mae","Data_points"];

    TT



