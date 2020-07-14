% Comments from annotations by Dan Welling (yellow page)
% Format of BATS-R-US output files
%
% Each value is padded by an integer (long, 32 bits) that lists how large
% the entry is in bits.
%
% The first entry is the ASCII header. It is padded as follows:
%
% [long int][header stuff in ascii][long int]
% ... where both long ints list the size of the header.
% The header contains units (not always?)
%
% The next entries (each line padded by length ints) are:
%
% - Simulation iteration (long int), simulation time (32 bit float), number
% of dimensions in file (long int), number of parameters (long int) and
% number of variables (long int).
%
% - Size of grid (number of points), a long for each dimension.
%
% - Parameters (a series of 32 bit floats).
%
% - ASCII list of variable + parameter names (a string of given length).
% The first values are the grid names (x, y, etc.), space separated. The
% next values are the variable names in order. The remaining values are the
% parameter names.
%
% - Vectors of the grid locations and variables, in order given by the
% previous entry.
%
% The first vector, usually 'x', is 32 bit floats x total number of points.

function [xmhd,ymhd,zmhd,vars] = read_BATSRUS_3dful(filename,type)
fid = fopen(filename);
%header_size = fread(fid,1,'long');
%header = fread(fid,header_size,'*char')';





%fid = fopen('shl_ful_4_n00050000.out');
% The first entry is the ASCII header. It is padded as follows:
%
% [long int][header stuff in ascii][long int]
% ... where both long ints list the size of the header.
% The header contains units (not always?)
header_size = fread(fid,1,'long');
header = fread(fid,header_size,'*char')';
header_size = fread(fid,1,'long');
% The next entries (each line padded by length ints) are:
%
% - Simulation iteration (long int), simulation time (32 bit float), number
% of dimensions in file (long int), number of parameters (long int) and
% number of variables (long int).
length1 = fread(fid,1,'long');
iter=fread(fid,1,'long');
simtime=fread(fid,1,'float');
numdim=fread(fid,1,'long');
numpar = fread(fid,1,'long');
numvar = fread(fid,1,'long');
length1 = fread(fid,1,'int');
% - Size of grid (number of points), a long for each dimension.
length1 = fread(fid,1,'int');
gridsizex = fread(fid,1,'long');
gridsizey = fread(fid,1,'long');
gridsizez = fread(fid,1,'long');
length1 = fread(fid,1,'int');
if(strcmp(type,'msmhd'))
    length1 = fread(fid,1,'int');
    param1 = fread(fid,1,'float');
    param2 = fread(fid,1,'float');
    param3 = fread(fid,1,'float');
    length1 = fread(fid,1,'int');
end
if(strcmp(type,'mfmhd'))
    length1 = fread(fid,1,'int');
    param1 = fread(fid,1,'float');
    param2 = fread(fid,1,'float');
    param3 = fread(fid,1,'float');
    param4 = fread(fid,1,'float');
    param5 = fread(fid,1,'float');
    param6 = fread(fid,1,'float');
    param7 = fread(fid,1,'float');
    length1 = fread(fid,1,'int');
end
% - ASCII list of variable + parameter names (a string of given length).
% The first values are the grid names (x, y, etc.), space separated. The
% next values are the variable names in order. The remaining values are the
% parameter names.
length1 = fread(fid,1,'int');
varspars = fread(fid,length1,'*char')';
varsparsaux=varspars;
disp(varspars)
disp(header)
count=size(varsparsaux,2);
while(varsparsaux(count) == ' ')
    varsparsaux(count) = [];
    count = size(varsparsaux,2);
end
blankspaces = strfind(varsparsaux,' ');
for count = 1:size(blankspaces,2);
    if count == 1
        varsparsaux(1:blankspaces(count)-1);
    elseif count == size(blankspaces,2)
        varsparsaux(blankspaces(count)+1:end);
    else
        varsparsaux(blankspaces(count-1)+1:blankspaces(count)-1);
    end
    
end
length1 = fread(fid,1,'int');
% - Vectors of the grid locations and variables, in order given by the
% previous entry.
npoints = gridsizex * gridsizey * gridsizez;
length1 = fread(fid,1,'int');
xmhd=fread(fid,npoints,'float');
ymhd=fread(fid,npoints,'float');
zmhd=fread(fid,npoints,'float');
length1 = fread(fid,1,'int');
if strcmp(type,'msmhd')
    dummy = 5;
elseif strcmp(type,'mfmhd')
    dummy = 9;
end
dummylongs = zeros(2,size(blankspaces,2)-dummy); % This should normally be 2, but it is five because the variable names include (xSI r g) that seems to be invalid.
vars = zeros(npoints,size(blankspaces,2)-dummy); % blankspaces + 1 (to capture all variables) - 3 (to account for the three already read)
for count = 1:size(blankspaces,2)-dummy
    dummylongs(1,count) = fread(fid,1,'int');
    vars(:,count) = fread(fid,npoints,'float');
    dummylongs(2,count) = fread(fid,1,'int');
end
fclose(fid);
end