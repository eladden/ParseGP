function [m,A,epoch,r,v] = generateList(filename, timefilecreated,maxdur, consts,whichconsts, minA, maxA)
%This function generates a list of satellites from an xml file.
%Usage:
% [m,A,epoch,r,v] = generateList(filename, consts)
% [m,A,epoch,r,v] = generateList(filename, consts, minA, maxA)
%
%Where:
%   filename - an xml file containing the TLE info
%   consts,whichconsts -  structures of constants for SGP4, can be generated with
%   generate_parameters function
%   minA, maxA - maximal and minimal area to generate random fake info
%
%   m - vector of length n containins the satellites' masses
%   A - a vector of length n containins the satellites' crossectional area 
%   (for drag and SRP)
%   epoch - times at which the state was taken
%   r - a matrix of nX3 containing the positions od the satellites
%   v - a matrix of nX3 containing the velocity of the satellites.
%

if nargin < 3
    error("Not enough input arguments")
elseif nargin < 4
    minA = 0.001;%m^2
    maxA = 0.1; %m^2
end

GPdata = readstruct(filename);

numberOfsats = length(GPdata.omm);

%preallocate memory
r = zeros(numberOfsats,3);
v = zeros(numberOfsats,3);
m = zeros(numberOfsats,1);
A = zeros(numberOfsats,1);
epoch = NaT(numberOfsats,1);

numOfDecayed = 0;
numOfStale = 0;
listCount = 1;

for i = 1:numberOfsats
    satstructxml = GPdata.omm(i).body.segment.data;
    satParamsNames = [satstructxml.userDefinedParameters.USER_DEFINED.parameterAttribute];
    decayPlace = find(satParamsNames == "DECAY_DATE");
    if isa(satstructxml.userDefinedParameters.USER_DEFINED(decayPlace).Text,"missing") 
        numOfDecayed = numOfDecayed +1;
        continue
    end
    if abs(timefilecreated - satstructxml.meanElements.EPOCH) > maxdur %The epoch was longer than maximal duration
        numOfStale = numOfStale + 1;
        continue
    end
    % Genetate the R and v
    satrec = GPxml2rv(whichconsts,consts,satstructxml);
    ep_ = satstructxml.meanElements.EPOCH;
    [~,r_,v_] = sgp4(satrec,0,consts);

    %check if the data has RCS_SIZE
    RCSPlace = find(satParamsNames == "RCS_SIZE");
    RCS_size = satstructxml.userDefinedParameters.USER_DEFINED(RCSPlace).Text;
    if isa(RCS_size,"string")
        switch satstructxml.userDefinedParameters.USER_DEFINED(RCSPlace).Text
            case "SMALL"
                minA = 0.01;
                maxA = 0.1;
            case "MEDIUM"
                minA = 0.1;
                maxA = 1;
            case "LARGE"
                minA = 1;
                maxA = 10;
        end
    end
    A_ = RCS2size(rand*(maxA-minA)+minA); %generate a random size
    %we want the A/m ratio to be maximum 1e-8 and minimum 1e-9 
    m_ = A_*1e-6/(rand*(1e-9-1e-8)+ 1e-8);
    
    %add to database
    r(listCount,:) = r_;
    v(listCount,:) = v_;
    A(listCount) = A_;
    m(listCount) = m_;
    epoch(listCount) = ep_;

    listCount = listCount+1;
end %end for
%
r = r(1:end-numOfDecayed-numOfStale,:);
v = v(1:end-numOfDecayed-numOfStale,:);
A = A(1:end-numOfDecayed-numOfStale)*1e-6; %km^2
m = m(1:end-numOfDecayed-numOfStale);
epoch = epoch(1:end-numOfDecayed-numOfStale);
end