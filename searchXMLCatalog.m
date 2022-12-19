function id = searchXMLCatalog(GPdata,value,searchtype,fieldname)
%id = searchXMLCatalog(GPdata,fieldname)
%search the structure containing the catalog of all satellites for value.
%You can search by field name by using fieldname variable, the default is 
% searching by NORAD caltalog ID
%You can search values that are greater than, less than, equal to or any of
%these kinds using 'eq','gt','lt','gte','lte' in the searchtype variable
%(the default is equal to)

if nargin < 4
    fieldname = 'NORAD_CAT_ID';
end
if nargin < 3
    searchtype = 'eq';
end

id = zeros(length(GPdata.omm));

switch fieldname
    case {'NORAD_CAT_ID','EPHEMERIS_TYPE','CLASSIFICATION_TYPE','ELEMENT_SET_NO'...
            'REV_AT_EPOCH','BSTAR','MEAN_MOTION_DOT','MEAN_MOTION_DDOT'}
        parentfield = 'tleParameters';
        grandparentfield = 'data';
    case {'EPOCH','MEAN_MOTION','ECCENTRICITY','INCLINATION','RA_OF_ASC_NODE',...
            'ARG_OF_PERICENTER','MEAN_ANOMALY'}
        parentfield = 'meanElements';
        grandparentfield = 'data';
    case {'OBJECT_NAME', 'OBJECT_ID', 'CENTER_NAME','REF_FRAME','TIME_SYSTEM',...
            'MEAN_ELEMENT_THEORY'}
        parentfield = 'metadata';
        if ~strcmp(searchtype,'eq')
            warning('The field you are searching is a text field, changing to equal search type');
        end
        searchtype = 'seq';
    otherwise
        error('unknown field')
end

for i = 1:length(GPdata.omm)
    switch searchtype
        case 'eq'
            if GPdata.omm(i).body.segment.(grandparentfield).(parentfield).(fieldname) == value
                id(i) = 1;
            end
        case 'lt'
            if GPdata.omm(i).body.segment.(grandparentfield).(parentfield).(fieldname) < value
                id(i) = 1;
            end
        case 'gt'
            if GPdata.omm(i).body.segment.(grandparentfield).(parentfield).(fieldname) > value
                id(i) = 1;
            end
        case 'gte'
            if GPdata.omm(i).body.segment.(grandparentfield).(parentfield).(fieldname) >= value
                id(i) = 1;
            end
        case 'lte'
            if GPdata.omm(i).body.segment.(parentfield).(fieldname) <= value
                id(i) = 1;
            end
        case 'eqs'
            if strcmp(GPdata.omm(i).body.segment.parentfield.fieldname,value)
                id(i) = 1;
            end
        otherwise
            error('huh?!')
    end
end

id = find(id);
end