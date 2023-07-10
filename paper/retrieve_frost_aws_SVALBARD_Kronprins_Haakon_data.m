

clear

clc
 opts = weboptions('ContentType', 'json','username','eddc3c74-7a9d-4c7a-a735-b16774458354');


 url      = 'https://frost.met.no/observations/v0.jsonld?';
 
 url_info = 'https://frost.met.no/sources/v0.jsonld?';
 



 varn = {'longitude','latitude','wind_speed','wind_from_direction','air_temperature','relative_humidity','surface_air_pressure'};
 vars = {'lon','lat','ws10','wd10','T2','RH2','P'};



% LIST OF SELECTED STATIONS BY STATION NUMBERS:


stnum = {'SN77035'}; % This is Kronprins Haakon


% https://frost.met.no/ex_userquest
stinfo = webread(url_info,opts,'types','SensorSystem','geometry','POLYGON ((5 74, 5 81, 40 81, 40 74))');
save stinfo stinfo


load stinfo



    validfrom = datestr([2021 03 01 00 00 00],'yyyy-mm-ddTHH');
    validto   = datestr([2022 11 01 00 00 00],'yyyy-mm-ddTHH');
    referencetime = [validfrom '/' validto];



AWS(1).name = 'Kronprins_Haakon';


st = 1;
    
%    try
    
    for vv = 1:length(varn)
    
        varnn = char(varn{vv});
        varss = char(vars{vv});
        
        varsst = [varss '_time'];


            rgb = webread(url,opts,'sources','sn77035','referencetime',referencetime,'elements',char(varn{vv}),'timeresolutions','PT10M');



                if vv == 1
%                   disp(num2str(st))  
                    disp(['Fetching ' AWS(st).name])
                end
                       
               for i = 1:length(rgb.data)
                AWS(st).(varsst)(i) = datenum(rgb.data(i).referenceTime,'yyyy-mm-ddTHH:MM:SS');
                AWS(st).(varss)(i) = rgb.data(i).observations(1).value;
               end
               
            
             
    end
        


for i = 1:length(AWS)

    DI = AWS(i).wd10;
    SP = AWS(i).ws10;

    AWS(i).u10 = SP.*cos((270-DI)*pi/180);
    AWS(i).v10 = SP.*sin((270-DI)*pi/180);

end

AWSN = AWS; clear AWS

a = intersect(AWSN(1).ws10_time,AWSN(1).lon_time);
b = intersect(AWSN(1).T2_time,a);
c = intersect(AWSN(1).RH2_time,b);
d = intersect(AWSN(1).P_time,c);



[a1 b1 c1] = intersect(AWSN(1).ws10_time,d);
    AWS(1).ws10 = AWSN(1).ws10(b1);
    AWS(1).wd10 = AWSN(1).wd10(b1);
   
[a1 b1 c1] = intersect(AWSN(1).T2_time,d);
    AWS(1).T2 = AWSN(1).T2(b1);

[a1 b1 c1] = intersect(AWSN(1).RH2_time,d);
    AWS(1).RH2 = AWSN(1).RH2(b1);

[a1 b1 c1] = intersect(AWSN(1).P_time,d);
    AWS(1).P = AWSN(1).P(b1);

[a1 b1 c1] = intersect(AWSN(1).lon_time,d);
    AWS(1).lon = AWSN(1).lon(b1);
    AWS(1).lat = AWSN(1).lat(b1);


AWS(1).time = AWSN(1).lon_time(b1);

AWS(1).mslp = AWS(1).P + 21/8;

    