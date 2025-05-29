function [m,A,epoch,r,v,satrecList] = generateList(filename, timefilecreated,maxdur, consts,whichconsts, minA, maxA)
%This function generates a list of satellites from an xml file.
%Usage:
% [m,A,epoch,r,v] = generateList(filename, timefilecreated,maxdur, consts,whichconsts)
% [m,A,epoch,r,v] = generateList(filename, timefilecreated,maxdur, consts,whichconsts, minA, maxA)
%
%Where:
%   filename - an xml file containing the TLE info
%   consts,whichconsts -  structures of constants for SGP4, can be generated with
%   generate_parameters function
%   timefilecreated - when was the xml file downloaded
%   minA, maxA - maximal and minimal area to generate random fake info
%
%   m - vector of length n containins the satellites' masses
%   A - a vector of length n containins the satellites' crossectional area 
%   (for drag and SRP)
%   epoch - times at which the state was taken
%   r - a matrix of nX3 containing the positions od the satellites
%   v - a matrix of nX3 containing the velocity of the satellites.
%   satrecList - a vector of structures (if you want to use SGP4 for
%   propagation)
%

if nargin < 3
    error("Not enough input arguments")
elseif nargin < 4
    minA = 0.001;%m^2
    maxA = 0.1; %m^2
end

tmpnames = [   {'error'        }
    {'satnum'       }
    {'jdsatepoch'   }
    {'ndot'         }
    {'nddot'        }
    {'bstar'        }
    {'inclo'        }
    {'nodeo'        }
    {'ecco'         }
    {'argpo'        }
    {'mo'           }
    {'no'           }
    {'a'            }
    {'alta'         }
    {'altp'         }
    {'isimp'        }
    {'method'       }
    {'aycof'        }
    {'con41'        }
    {'cc1'          }
    {'cc4'          }
    {'cc5'          }
    {'d2'           }
    {'d3'           }
    {'d4'           }
    {'delmo'        }
    {'eta'          }
    {'argpdot'      }
    {'omgcof'       }
    {'sinmao'       }
    {'t'            }
    {'t2cof'        }
    {'t3cof'        }
    {'t4cof'        }
    {'t5cof'        }
    {'x1mth2'       }
    {'x7thm1'       }
    {'mdot'         }
    {'nodedot'      }
    {'xlcof'        }
    {'xmcof'        }
    {'nodecf'       }
    {'irez'         }
    {'d2201'        }
    {'d2211'        }
    {'d3210'        }
    {'d3222'        }
    {'d4410'        }
    {'d4422'        }
    {'d5220'        }
    {'d5232'        }
    {'d5421'        }
    {'d5433'        }
    {'dedt'         }
    {'del1'         }
    {'del2'         }
    {'del3'         }
    {'didt'         }
    {'dmdt'         }
    {'dnodt'        }
    {'domdt'        }
    {'e3'           }
    {'ee2'          }
    {'peo'          }
    {'pgho'         }
    {'pho'          }
    {'pinco'        }
    {'plo'          }
    {'se2'          }
    {'se3'          }
    {'sgh2'         }
    {'sgh3'         }
    {'sgh4'         }
    {'sh2'          }
    {'sh3'          }
    {'si2'          }
    {'si3'          }
    {'sl2'          }
    {'sl3'          }
    {'sl4'          }
    {'gsto'         }
    {'xfact'        }
    {'xgh2'         }
    {'xgh3'         }
    {'xgh4'         }
    {'xh2'          }
    {'xh3'          }
    {'xi2'          }
    {'xi3'          }
    {'xl2'          }
    {'xl3'          }
    {'xl4'          }
    {'xlamo'        }
    {'zmol'         }
    {'zmos'         }
    {'atime'        }
    {'xli'          }
    {'xni'          }
    {'operationmode'}
    {'no_kozai'     }
    {'am'           }
    {'em'           }
    {'im'           }
    {'Om'           }
    {'mm'           }
    {'nm'           }
    {'tumin'        }
    {'mu'           }
    {'radiusearthkm'}
    {'xke'          }
    {'j2'           }
    {'j3'           }
    {'j4'           }
    {'j3oj2'        }
    {'init'         }];
tmpvalues = repmat({0},115,1);
tmp = cell2struct(tmpvalues,tmpnames);

GPdata = readstruct(filename);

numberOfsats = length(GPdata.omm);

%preallocate memory
r = zeros(numberOfsats,3);
v = zeros(numberOfsats,3);
m = zeros(numberOfsats,1);
A = zeros(numberOfsats,1);
epoch = NaT(numberOfsats,1);
satrecList = repmat(tmp,numberOfsats,1);

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
    ep_ = satstructxml.meanElements.EPOCH;
    [~,r_,v_] = sgp4(satrec,0);%,consts);

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
    satrecList(listCount) = satrec;
    listCount = listCount+1;
end %end for
%
r = r(1:end-numOfDecayed-numOfStale,:);
v = v(1:end-numOfDecayed-numOfStale,:);
A = A(1:end-numOfDecayed-numOfStale)*1e-6; %km^2
m = m(1:end-numOfDecayed-numOfStale);
epoch = epoch(1:end-numOfDecayed-numOfStale);
satrecList = satrecList(1:end-numOfDecayed-numOfStale);
end