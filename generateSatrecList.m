function [satrecList] = generateSatrecList(filename, timefilecreated,maxdur, consts,whichconsts)
%This function generates a list of satellites from an xml file.
%Usage:
% [satrecList] = generateSatrecList(filename, timefilecreated,maxdur, consts,whichconsts, minA, maxA)
%
%Where:
%   filename - an xml file containing the TLE info
%   consts,whichconsts -  structures of constants for SGP4, can be generated with
%   generate_parameters function
%   minA, maxA - maximal and minimal area to generate random fake info
%
%   satrecList - a vector of structures (if you want to use SGP4 for
%   propagation)
%

GPdata = readstruct(filename);

numberOfsats = length(GPdata.omm);

%preallocate memory

%satrecList = repmat(struct('a',0),numberOfsats,1);

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
    satrec = GPxml2rv(whichconsts,'i',consts,satstructxml);
    satrecList(listCount) = satrec;
    listCount = listCount+1;
end %end for
%
%satrecList = satrecList(1:end-numOfDecayed-numOfStale);
end