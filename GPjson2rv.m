function [satrec, startmfe, stopmfe, deltamin] = GPxml2rv(whichconst, opsmode,consts, struct)
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
    xpdotp   =  1440.0 / (2.0*pi);   %  229.1831180523293;  % [rev/day]/[rad/min]  

 
    satrec.error = 0;

    satrec.satnum = struct.NORAD_CAT_ID;
    satrec.jdsatepoch = juliandate(datetime(struct.EPOCH));
    satrec.ndot = struct.MEAN_MOTION_DOT / (xpdotp*1440.0); % [rad/min^2]
    satrec.nddot = struct.MEAN_MOTION_DDOT / (xpdotp*1440.0*1440); % [rad/min^3]
    satrec.bstar = struct.BSTAR;
    satrec.inclo = deg2rad(struct.INCLINATION);
    satrec.nodeo = deg2rad(struct.RA_OF_ASC_NODE);
    satrec.ecco = struct.ECCENTRICITY;
    satrec.argpo = deg2rad(struct.ARG_OF_PERICENTER);
    satrec.mo = deg2rad(struct.MEAN_ANOMALY);
    satrec.no =  struct.MEAN_MOTION / xpdotp;

    tumin = consts(1);

    satrec.a    = (satrec.no*tumin)^(-2/3);

    satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;



    sgp4epoch = satrec.jdsatepoch - 2433281.5; % days since 0 Jan 1950
    [satrec] = sgp4init( whichconst, opsmode, satrec, satrec.jdsatepoch-2433281.5, satrec.bstar, ...
        satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no, ...
        satrec.nodeo);