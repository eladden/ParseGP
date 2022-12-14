function [satrec, startmfe, stopmfe, deltamin] = GPxml2rv(whichconst, tumin, struct)
%This function replaces the old twoline2rv.m that was used to read two
%strings of TLE and produced the satrec structure used the the SGP4
%implementation.
%to use this function you need to read an xml data file with
%
%GPdata=readstruct(filename) then use satGP = GPdata.
%
%Then access the satellite using
%
%sat_num_i_struct = GPxml.omm(i).body.segment.data

    %global tumin;

    %deg2rad  =   pi / 180.0;         %  0.01745329251994330;  % [deg/rad]
    xpdotp   =  1440.0 / (2.0*pi);   % 229.1831180523293;  % [rev/day]/[rad/min]  

 
    satrec.error = 0;

    satrec.satnum = struct.tleParameters.NORAD_CAT_ID;
    satrec.jdsatepoch = juliandate(struct.meanElements.EPOCH);
    satrec.ndot = struct.tleParameters.MEAN_MOTION_DOT;
    satrec.nddot = struct.tleParameters.MEAN_MOTION_DDOT;
    satrec.bstar = struct.tleParameters.BSTAR;
    satrec.inclo = deg2rad(struct.meanElements.INCLINATION);
    satrec.nodeo = deg2rad(struct.meanElements.RA_OF_ASC_NODE);
    satrec.ecco = struct.meanElements.ECCENTRICITY;
    satrec.argpo = deg2rad(struct.meanElements.ARG_OF_PERICENTER);
    satrec.mo = deg2rad(struct.meanElements.MEAN_ANOMALY);
    satrec.no =  struct.meanElements.MEAN_MOTION * xpdotp;

    satrec.a    = (satrec.no*tumin)^(-2/3);

    satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;



    sgp4epoch = satrec.jdsatepoch - 2433281.5; % days since 0 Jan 1950
    [satrec] = sgp4init(whichconst, satrec, satrec.bstar, satrec.ecco, sgp4epoch, ...
         satrec.argpo, satrec.inclo, satrec.mo, satrec.no, satrec.nodeo);