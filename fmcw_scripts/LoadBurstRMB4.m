function vdat = LoadBurstRMB4(Filename, Hdr, Burst, SamplingFrequency)

% vdat = LoadBurstRMB4(Filename, Burst, SamplesPerChirp)
%
% Read FMCW data file from RMB-B? (Data from after Oct 2013)

% 2014/5/22 Add temperature correction CLS

% Updated by Falk Oraschewski, 20.03.2021, changeg strfind()
% find(ismember()), because strfind() gave an empty return.

WperChirpHdr = 0;
MaxHeaderLen = length(Hdr);
burstpointer = 0;
vdat.Code = 0;
fid = fopen(Filename,'r');
if fid >= 0
    fseek(fid,0,'eof');
    filelength = ftell(fid);
    BurstCount = 1;
    while BurstCount <= Burst && burstpointer <= filelength - MaxHeaderLen
        fseek(fid,burstpointer,'bof');
        searchind = strfind(Hdr, 'Samples:');
        if ~isempty(searchind)
            try
                searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
                vdat.Nsamples = sscanf(Hdr(searchind(1)+8:searchCR(1)+searchind(1)),'%d');
                WperChirpCycle = vdat.Nsamples + WperChirpHdr;
                searchind = strfind(Hdr, 'SubBursts in burst:');
                searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
                vdat.SubBurstsInBurst = sscanf(Hdr(searchind(1)+19:searchCR(1)+searchind(1)),'%d');

                searchind = strfind(Hdr, 'Average:');
                if isempty(searchind)
                    vdat.Average = 0; %cls 9/jan/14 -average not included in mooring deploy
                else
                searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
                    vdat.Average = sscanf(Hdr(searchind(1)+8:searchCR(1)+searchind(1)),'%d');
                end
                
                searchind = strfind(Hdr, 'nAttenuators:');
                searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
                vdat.NAttenuators = sscanf(Hdr(searchind(1)+13:searchCR(1)+searchind(1)),'%d',1);

                searchind = strfind(Hdr, 'Attenuator 1:');
                searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
                vdat.Attenuator_1 = sscanf(Hdr(searchind(1)+13:searchCR(1)+searchind(1)),'%f,',vdat.NAttenuators);
                
                searchind = strfind(Hdr, 'Attenuator 2:');
                searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
                vdat.Attenuator_2 = sscanf(Hdr(searchind(1)+13:searchCR(1)+searchind(1)),'%f,',vdat.NAttenuators);

                searchind = strfind(Hdr, 'Tx Antenna select:');
                searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
                vdat.TxAnt = sscanf(Hdr(searchind(1)+18:searchCR(1)+searchind(1)),'%d',8);
                
                searchind = strfind(Hdr, 'Rx Antenna select:');
                searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
                vdat.RxAnt = sscanf(Hdr(searchind(1)+18:searchCR(1)+searchind(1)),'%d',8);
                
                ind = find(vdat.TxAnt~=1);
                vdat.TxAnt(ind) = [];
                ind = find(vdat.RxAnt~=1);
                vdat.RxAnt(ind) = [];
                
                if vdat.Average % commented out by Craig Stewart 18 March 2014 as this causes problems if running in attended mode with average set to 1 (all chirps are present but were lost due to setting ChirpsInBurst to 1)
                    vdat.ChirpsInBurst = 1;
                else
                    vdat.ChirpsInBurst = vdat.SubBurstsInBurst * length(vdat.TxAnt) * ...
                       length(vdat.RxAnt) * vdat.NAttenuators;
               end
                               
                searchind = strfind(Hdr, '*** End Header ***');
                
                burstpointer = burstpointer + searchind(1) + 20;
            catch
                vdat.Code = -2;
                vdat.Burst = BurstCount;
                %keyboard
                return
            end
        end
        WordsPerBurst = vdat.ChirpsInBurst * WperChirpCycle;
        if BurstCount < Burst && burstpointer <= filelength - MaxHeaderLen
            if vdat.Average == 2
                burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle*4;
            else
                burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle*2;
            end
        end
        BurstCount = BurstCount + 1;
    end
    
    % Extract remaining information from header
    searchind = strfind(Hdr, 'Time stamp:');
    if isempty(searchind)
        vdat.Code = -4;
        return
    end
    try
        searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
        td = sscanf(Hdr(searchind(1)+11:searchCR(1)+searchind(1)),...
            '%d-%d-%d %d:%d:%d');
        vdat.TimeStamp = datenum(td(1),td(2),td(3),td(4),td(5),td(6));
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(Hdr, 'Temperature 1:');
    try
        searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
        vdat.Temperature_1 = sscanf(Hdr(searchind(1)+14:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(Hdr, 'Temperature 2:');
    try
        searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
        vdat.Temperature_2 = sscanf(Hdr(searchind(1)+14:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(Hdr, 'Battery voltage:');
    try
        searchCR = find(ismember(Hdr(searchind(1):end), char([10 13])));
        vdat.BatteryVoltage = sscanf(Hdr(searchind(1)+16:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    
    fseek(fid,burstpointer-1,'bof');
    if BurstCount == Burst+1
        if vdat.Average == 2
            [vdat.v count] = fread(fid,WordsPerBurst,'*uint32','ieee-le');
        else
            [vdat.v count] = fread(fid,WordsPerBurst,'*uint16','ieee-le');
        end
        if count < WordsPerBurst
            vdat.Code = 2;
        end
        vdat.v(vdat.v<0) = vdat.v(vdat.v<0) + 2^16;
        vdat.v = single(vdat.v);
        vdat.v = vdat.v * 2.5 / 2^16;
        if vdat.Average == 2
            vdat.v = vdat.v / (vdat.SubBurstsInBurst * vdat.NAttenuators);
        end
        vdat.Startind = ((WperChirpHdr+1):WperChirpCycle:WperChirpCycle*vdat.ChirpsInBurst)';
        vdat.Endind = vdat.Startind + WperChirpCycle - 1;
        vdat.Burst = Burst;
    else
        % Too few bursts in file
        vdat.Burst = BurstCount - 1;
        vdat.Code = -4;
        %keyboard
    end
    fclose(fid);
else
    % Unknown file
    vdat.Code = -1;
end

% Clean temperature record (wrong data type?)
bti1 = find(vdat.Temperature_1>300); % bad temperature indices
if ~isempty(bti1)
    %disp('Cleaning temperature over 300C')
    vdat.Temperature_1(bti1) = vdat.Temperature_1(bti1)-512;
end
bti2 = find(vdat.Temperature_2>300); % bad temperature indices
vdat.Temperature_2(bti2) = vdat.Temperature_2(bti2)-512;
