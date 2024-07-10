function [fmt,Hdr] = fmcw_file_format(FilePath)

% fmt = fmcw_file_format(filename)
%
% Determine fmcw file format from burst header using keyword presence

% Craig Stewart
% 2013-10-20
%
% Updated by Keith Nicholls, 2014-10-22: RMB2
%
% Updated by Reza Ershadi, Feb 2021: use fscanf to read and extract the
% header. fread had some problem with some the ApRES data.

[fid, msg] = fopen(FilePath,'r');
if fid == -1
    error(msg)
end
A = fscanf(fid,'%c');
%A1 = fread(fid,2000,'*char');
fclose(fid);
HeaderEndFlag = '*** End Header ***';
a1 = strfind(A,HeaderEndFlag);
MaxHeaderLen = a1+length(HeaderEndFlag);
Hdr = A(1:MaxHeaderLen);

if contains(Hdr, 'IQ=1', 'IgnoreCase',true) % IQ data
    fmt = 6;
elseif contains(Hdr, 'SW_Issue=', 'IgnoreCase',true) % Data from RMB2 after Oct 2014
    fmt = 5;
elseif contains(Hdr, 'SubBursts in burst:', 'IgnoreCase',true) % Data from after Oct 2013
    fmt = 4;
elseif contains(Hdr, '*** Burst Header ***', 'IgnoreCase',true) % Data from Jan 2013
    fmt = 3;
elseif contains(Hdr, 'RADAR TIME', 'IgnoreCase',true) % Data from Prototype FMCW radar (nov 2012)
    fmt = 2;
else
    %fmt = 0; % unknown file format
    error('Unknown file format - check file')
end
