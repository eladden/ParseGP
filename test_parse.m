%Test the parsing function
%% general parameters
fprintf('loading...');
if ~exist('consts','var')
    [xpdotp, Re, J2, mu, whichconsts, consts, finals, dAT] = generate_parameters;
end
tumin = consts(1);
fprintf('done\n');
%% generate the old way - using TLE
%SOYUZ-MS 27             
str1 = '1 63520U 25072A   25127.51095612  .00007625  00000+0  14445-3 0  9998';
str2 = '2 63520  51.6346 148.5519 0002315  86.8352 273.2901 15.49392895  2482';

[~,~,~,satrec_old] =  twoline2rv(str1, str2, 'u', 'm','i',whichconsts);%


%% generate using GP data from xml
GPdata = readstruct("last-30-days.xml");
satstructxml = GPdata.omm(1).body.segment.data;

satrec_new_xml = GPxml2rv(whichconsts,'i',consts,satstructxml);

fn = fieldnames(satrec_new_xml);
for i = 1 : numel(fn)
  errxml.(fn{i}) = abs(satrec_new_xml.(fn{i}) - satrec_old.(fn{i}));
end
%%
GPdata = jsondecode(fileread("last-30-days.json"));
satstructjson = GPdata(1);

satrec_new_xml = GPjson2rv(whichconsts,'i',consts,satstructjson);

for i = 1 : numel(fn)
  errjson.(fn{i}) = abs(satrec_new_xml.(fn{i}) - satrec_old.(fn{i}));
end
%