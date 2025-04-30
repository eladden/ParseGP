%Test the parsing function
%% general parameters
fprintf('loading...');
if ~exist('consts','var')
    [xpdotp, Re, J2, mu, whichconsts, consts, finals, dAT] = generate_parameters;
end
tumin = consts(1);
fprintf('done\n');
%% generate the old way - using TLE
%YAOGAN-34 03    
str1 = '1 54249U 22154A   22348.16127061 -.00000090  00000+0  00000+0 0  9999';
str2 = '2 54249  63.4088 311.5516 0008941 272.2943  87.7052 13.45255822  3928';

satrec_old =  twoline2rv(whichconsts, str1, str2, 'c', [], consts);%


%% generate using GP data from xml
GPdata = readstruct("last-30-days.xml");
satstructxml = GPdata.omm(1).body.segment.data;

satrec_new = GPxml2rv(whichconsts,consts,satstructxml);

fn = fieldnames(satrec_new);
for i = 1 : numel(fn)
  err.(fn{i}) = abs(satrec_new.(fn{i}) - satrec_old.(fn{i}));
end

%