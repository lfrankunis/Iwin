%% Matching boats with Narveneset


IWIN(4).lat = 78.56343;
IWIN(4).lon = 16.29687;

for id = 1:3
    id
    for i = 1:length(IWIN(id).lon)
           IWIN(id).NNdist(i) = m_lldist([IWIN(4).lon IWIN(id).lon(i)],[IWIN(4).lat IWIN(id).lat(i)]);
    end
end



for id = 1:3


    a = zeros(length(IWIN(id).lon),1); b=a; c=a; d=a;


    ida = datenum([2020 06 18 00 00 00]);
    idb = datenum([2023 06 23 00 00 00]);

    a(IWIN(id).time >= ida & IWIN(id).time <= idb)  = 1;
    b(IWIN(id).GPS_speed > 0.25) = 1;
    c(IWIN(id).NNdist < 2) = 1;

dr = a + b + c;

tt = find(dr == 3);


   IWIN(id).temperature(IWIN(id).exhaust_plume_influence == 1) = nan;


            iwintime =  IWIN(id).time(tt);



       [a b c] = intersect(iwintime,IWIN(4).time);


        IWINs(id).time = IWIN(id).time(tt(b));
        IWINs(id).lon = IWIN(id).lon(tt(b));
        IWINs(id).lat = IWIN(id).lat(tt(b));
        IWINs(id).temperature = IWIN(id).temperature(tt(b));
        IWINs(id).relative_humidity = IWIN(id).relative_humidity(tt(b));
        IWINs(id).air_pressure = IWIN(id).air_pressure(tt(b));
        IWINs(id).mslp = IWIN(id).mslp(tt(b));
        IWINs(id).wind_speed_corrected = IWIN(id).wind_speed_corrected(tt(b));
        IWINs(id).wind_direction_corrected = IWIN(id).wind_direction_corrected(tt(b));
        IWINs(id).GPS_heading = IWIN(id).GPS_heading(tt(b));
        IWINs(id).GPS_speed = IWIN(id).GPS_speed(tt(b));

        IWINs_NN(id).time = IWIN(4).time(c);
        IWINs_NN(id).temperature = IWIN(4).temperature(c);
        IWINs_NN(id).relative_humidity = IWIN(4).relative_humidity(c);
        IWINs_NN(id).air_pressure = IWIN(4).air_pressure(c);
        IWINs_NN(id).mslp = IWIN(4).mslp(c);
        IWINs_NN(id).wind_speed = IWIN(4).wind_speed(c);
        IWINs_NN(id).wind_direction = IWIN(4).wind_direction(c);


end


     %%

    clc
    clear *_diff *_mae *_points
    for i = 1:3
        temp_diff(i) = nanmean(IWINs(i).temperature - IWINs_NN(i).temperature);
        temp_mae(i) = nanmean(abs(IWINs(i).temperature - IWINs_NN(i).temperature));

        rh_diff(i) = nanmean(IWINs(i).relative_humidity - IWINs_NN(i).relative_humidity);
        rh_mae(i) = nanmean(abs(IWINs(i).relative_humidity - IWINs_NN(i).relative_humidity));

        p_diff(i) = nanmean(IWINs(i).air_pressure - IWINs_NN(i).air_pressure);
        p_mae(i) = nanmean(abs(IWINs(i).air_pressure - IWINs_NN(i).air_pressure));

        mslp_diff(i) = nanmean(IWINs(i).mslp - IWINs_NN(i).mslp);
        mslp_mae(i) = nanmean(abs(IWINs(i).mslp - IWINs_NN(i).mslp));

        ws_diff(i) = nanmean(IWINs(i).wind_speed_corrected - IWINs_NN(i).wind_speed);
        ws_mae(i) = nanmean(abs(IWINs(i).wind_speed_corrected - IWINs_NN(i).wind_speed));

            a11 = IWINs(i).wind_direction_corrected;
            b11 = IWINs_NN(i).wind_direction;
            a = rad2deg(angdiff(deg2rad(a11),deg2rad(b11)));
            a = a(isfinite(a));
            a_abs = abs(a);
        wd_diff(i) = rad2deg(circ_mean(deg2rad(a),[],1));
        wd_mae(i) = rad2deg(circ_mean(deg2rad(a_abs),[],1));

        data_points(i) = length(IWINs(i).temperature);
    end


    temp_diffs = round(temp_diff,2);
    rh_diffs = round(rh_diff,2);
    p_diffs = round(p_diff,2);
    mslp_diffs = round(mslp_diff,2);
    ws_diffs = round(ws_diff,2);
    wd_diffs = round(wd_diff,2);

    temp_maes = round(temp_mae,2);
    rh_maes = round(rh_mae,2);
    p_maes = round(p_mae,2);
    mslp_maes = round(mslp_mae,2);
    ws_maes = round(ws_mae,2);
    wd_maes = round(wd_mae,2);


%%

clc
clear T
Comparisons = ["Bard vs Narveneset","Polargirl vs Narveneset","Billefjord vs Narveneset"];

  T = table(Comparisons',temp_diffs',rh_diffs',mslp_diffs',ws_diffs',wd_diffs',data_points');
    T.Properties.VariableNames = ["Comparison","Temperature_diff","RelHumidity_diff","Mslp_diff","Windspeed_diff","Winddirection_diff","Data_points"];

    T


%%

clc
clear TT
Comparisons = ["Bard vs Narveneset","Polargirl vs Narveneset","Billefjord vs Narveneset"];

  TT = table(Comparisons',temp_maes',rh_maes',mslp_maes',ws_maes',wd_maes',data_points');
    TT.Properties.VariableNames = ["Comparison","Temperature_mae","RelHumidity_mae","Mslp_mae","Windspeed_mae","Winddirection_mae","Data_points"];

    TT



