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
    satrecList(listCount) = satrec;
    listCount = listCount+1;
end %end for
%
satrecList = satrecList(1:end-numOfDecayed-numOfStale);
end