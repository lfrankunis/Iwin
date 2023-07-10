

clear

clc
 opts = weboptions('ContentType', 'json','username','eddc3c74-7a9d-4c7a-a735-b16774458354');

 url      = 'https://frost.met.no/observations/v0.jsonld?';
 
 url_info = 'https://frost.met.no/sources/v0.jsonld?';
 


 varn = {'wind_speed','wind_from_direction','air_temperature','relative_humidity','surface_air_pressure'};
 vars = {'ws10','wd10','T2','RH2','P'};



% LIST OF SELECTED STATIONS BY STATION NUMBERS:
stnums = {'SN99840'}; % This is Svalbard airport


% https://frost.met.no/ex_userquest
stinfo = webread(url_info,opts,'types','SensorSystem','geometry','POLYGON ((5 74, 5 81, 40 81, 40 74))');
save stinfo stinfo


load stinfo


for i = 1:length(stinfo.data)
   stname{i} =  stinfo.data{i,1}.id;
end
for i = 1:length(stnums)
    stn=stnums{i};
    stnum{i}=str2num(stn(3:end));
    
    kk = find(strcmp(stnums{i}, stname));
    AWS(i).name = stinfo.data{kk,1}.name;
    AWS(i).lon  = stinfo.data{kk,1}.geometry.coordinates(1);
    AWS(i).lat  = stinfo.data{kk,1}.geometry.coordinates(2);
end



yrs = 2021:2023;


for st = 1:length(stnum)
    
    
    for vv = 1:length(varn)
    
        varnn = char(varn{vv});
        varss = char(vars{vv});
        
        varsst = [varss '_time'];


        AWS(st).(varsst) = [];
        AWS(st).(varss) = [];

        for qq = 1:length(yrs)

                validfrom = datestr([yrs(qq) 01 01 00 00 00],'yyyy-mm-ddTHH');
                validto   = datestr([yrs(qq) 12 31 00 00 00],'yyyy-mm-ddTHH');
                referencetime = [validfrom '/' validto];


                %         try 
                if strcmp(varss,'ws10') || strcmp(varss,'wd10')
                            rgb = webread(url,opts,'sources',['sn' num2str(stnum{st})],'referencetime',referencetime,'elements',char(varn{vv}),'timeresolutions','PT10M');
                elseif strcmp(varss,'T2') || strcmp(varss,'RH2') || strcmp(varss,'P')
                            rgb = webread(url,opts,'sources',['sn' num2str(stnum{st})],'referencetime',referencetime,'elements',char(varn{vv}),'timeresolutions','PT1H');
                end


                if vv == 1
                    disp(num2str(st))  
                    disp(['Fetching ' AWS(st).name])
                end
                       
               
 

               clear aa bb
               for i = 1:length(rgb.data)
                aa(i) = datenum(rgb.data(i).referenceTime,'yyyy-mm-ddTHH:MM:SS');
                bb(i) = rgb.data(i).observations(1).value;
               end

        AWS(st).(varsst) = [AWS(st).(varsst) aa];
        AWS(st).(varss)  = [AWS(st).(varss) bb];


        end

       
    end
        


end


AWSN = AWS; clear AWS

a = intersect(AWSN(1).ws10_time,AWSN(1).wd10_time);
[a1 b1 c1] = intersect(AWSN(1).ws10_time,a);
    AWS(1).ws10 = AWSN(1).ws10(b1);

[a1 b1 c1] = intersect(AWSN(1).wd10_time,a);
    AWS(1).wd10 = AWSN(1).wd10(b1);

    AWS(1).wind_time = AWSN(1).wd10_time(b1);


b = intersect(AWSN(1).T2_time,AWSN(1).RH2_time);
c = intersect(AWSN(1).P_time,b);

[a1 b1 c1] = intersect(AWSN(1).T2_time,c);
    AWS(1).T2 = AWSN(1).T2(b1);

[a1 b1 c1] = intersect(AWSN(1).RH2_time,c);
    AWS(1).RH2 = AWSN(1).RH2(b1);

[a1 b1 c1] = intersect(AWSN(1).P_time,c);
    AWS(1).P = AWSN(1).P(b1);

AWS(1).PTH_time = AWSN(1).P_time(b1);

AWS(1).lon = AWSN(1).lon;
AWS(1).lat = AWSN(1).lat;



AWS(1).mslp = AWS(1).P + 28/8;

