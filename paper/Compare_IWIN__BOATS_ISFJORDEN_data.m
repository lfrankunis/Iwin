



%% Matching boat locations on Isfjorden



clear Temperature_diff Humidity_diff Windspeed_diff Winddirection_diff

for ii = 1:3


    if ii == 1
        idf = [1 2];
    elseif ii == 2
        idf = [1 3];
    elseif ii == 3
        idf = [2 3];
    end


ida = datenum([2021 03 01 00 00 00]);
idb = datenum([2023 06 23 00 00 00]);

clc
for id = idf

    a = zeros(length(IWIN(id).lon),1); b=a; c=a; d=a;


        a(IWIN(id).lon < 20 & IWIN(id).lon > 11.9) = 1;
        b(IWIN(id).lat > 77 & IWIN(id).lat < 79) = 1;
        c(IWIN(id).GPS_speed > 0.25) = 1;


    IWIN(id).temperature(IWIN(id).exhaust_plume_influence == 1) = nan; 


    dd = a + b + c;


    ids(id).n = find(dd == 3);
end


k1 = intersect(IWIN(idf(1)).time(ids(idf(1)).n),IWIN(idf(2)).time(ids(idf(2)).n));

for id = idf
    [a b c] = intersect(k1,IWIN(id).time);
    IWINs(id).time = IWIN(id).time(c);
    IWINs(id).lon = IWIN(id).lon(c);
    IWINs(id).lat = IWIN(id).lat(c);
    IWINs(id).temperature = IWIN(id).temperature(c);
    IWINs(id).relative_humidity = IWIN(id).relative_humidity(c);
    IWINs(id).mslp = IWIN(id).mslp(c);
    IWINs(id).wind_speed_corrected = IWIN(id).wind_speed_corrected(c);
    IWINs(id).wind_direction_corrected = IWIN(id).wind_direction_corrected(c);
    IWINs(id).GPS_speed = IWIN(id).GPS_speed(c);
    IWINs(id).GPS_heading = IWIN(id).GPS_heading(c);
end

clear disst
for i = 1:length(IWINs(idf(1)).time)
    disst(i) = m_lldist([IWINs(idf(1)).lon(i) IWINs(idf(2)).lon(i)],[IWINs(idf(1)).lat(i) IWINs(idf(2)).lat(i)]);
end


  kk = find(disst < 1);

  mp(ii).kk = kk;
  mp(ii).idf = idf;
  mp(ii).lon_1 = IWINs(idf(1)).lon(kk); 
  mp(ii).lon_2 = IWINs(idf(2)).lon(kk);
  mp(ii).lat_1 = IWINs(idf(1)).lat(kk); 
  mp(ii).lat_2 = IWINs(idf(2)).lat(kk);

  mp(ii).temperature_1 = IWINs(idf(1)).temperature(kk);
  mp(ii).temperature_2 = IWINs(idf(2)).temperature(kk);
  mp(ii).relative_humidity_1 = IWINs(idf(1)).relative_humidity(kk);
  mp(ii).relative_humidity_2 = IWINs(idf(2)).relative_humidity(kk);
  mp(ii).mslp_1 = IWINs(idf(1)).mslp(kk);
  mp(ii).mslp_2 = IWINs(idf(2)).mslp(kk);
  mp(ii).wind_speed_corrected_1 = IWINs(idf(1)).wind_speed_corrected(kk);
  mp(ii).wind_speed_corrected_2 = IWINs(idf(2)).wind_speed_corrected(kk);
  mp(ii).wind_direction_corrected_1 = IWINs(idf(1)).wind_direction_corrected(kk);
  mp(ii).wind_direction_corrected_2 = IWINs(idf(2)).wind_direction_corrected(kk);

  mp(ii).name_1 = IWIN(idf(1)).name;
  mp(ii).name_2 = IWIN(idf(2)).name;

    temp_diff(ii) = nanmean(IWINs(idf(1)).temperature(kk) - IWINs(idf(2)).temperature(kk));
    rh_diff(ii) = nanmean(IWINs(idf(1)).relative_humidity(kk) - IWINs(idf(2)).relative_humidity(kk));
    mslp_diff(ii) = nanmean(IWINs(idf(1)).mslp(kk) - IWINs(idf(2)).mslp(kk));
    ws_diff(ii) = nanmean(IWINs(idf(1)).wind_speed_corrected(kk) - IWINs(idf(2)).wind_speed_corrected(kk));


    temp_mae(ii) = nanmean(abs(IWINs(idf(1)).temperature(kk) - IWINs(idf(2)).temperature(kk)));
    rh_mae(ii) = nanmean(abs(IWINs(idf(1)).relative_humidity(kk) - IWINs(idf(2)).relative_humidity(kk)));
    mslp_mae(ii) = nanmean(abs(IWINs(idf(1)).mslp(kk) - IWINs(idf(2)).mslp(kk)));
    ws_mae(ii) = nanmean(abs(IWINs(idf(1)).wind_speed_corrected(kk) - IWINs(idf(2)).wind_speed_corrected(kk)));


a11 = IWINs(idf(1)).wind_direction_corrected(kk);
b11 = IWINs(idf(2)).wind_direction_corrected(kk);

        a = rad2deg(angdiff(deg2rad(a11),deg2rad(b11)));
        a = a(isfinite(a));
            a_abs = abs(a);
        wd_diff(ii) = rad2deg(circ_mean(deg2rad(a),[],1));
        wd_mae(ii) = rad2deg(circ_mean(deg2rad(a_abs),[],1));


    Data_points(ii) = length(kk);

end


    temp_diffs = round(temp_diff,2);
    rh_diffs   = round(rh_diff,2);
%     p_diffs    = round(p_diff,2);
    mslp_diffs = round(mslp_diff,2);
    ws_diffs   = round(ws_diff,2);
    wd_diffs   = round(wd_diff,2);

    temp_maes = round(temp_mae,2);
    rh_maes   = round(rh_mae,2);
%     p_maes    = round(p_mae,2);
    mslp_maes = round(mslp_mae,2);
    ws_maes   = round(ws_mae,2);
    wd_maes   = round(wd_mae,2);

%%
clc
clear T
Comparisons = ["Bard vs Polargirl","Bard vs Billefjord","Polargirl vs Billefjord"];

  T = table(Comparisons',temp_diffs',rh_diffs',mslp_diffs',ws_diffs',wd_diffs',Data_points');
    T.Properties.VariableNames = ["Comparison","Temperature_diff","RelHumidity_diff","Mslp_diff","Windspeed_diff","Winddirection_diff","Data_points"];

    T

%%
clc
    clear TT
Comparisons = ["Bard vs Polargirl","Bard vs Billefjord","Polargirl vs Billefjord"];

  TT = table(Comparisons',temp_maes',rh_maes',mslp_maes',ws_maes',wd_maes',Data_points');
    TT.Properties.VariableNames = ["Comparison","Temperature_mae","RelHumidity_mae","Mslp_mae","Windspeed_mae","Winddirection_mae","Data_points"];

    TT





