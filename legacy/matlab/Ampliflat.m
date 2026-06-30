function ampliflat(x,atype,options)

%AMPLIFLAT Ideal Optical amplifier with ASE noise.
%   AMPLIFLAT(X,ATYPE) amplifies the optical field. ATYPE is a string     
%   equal to 'gain' if the amplifier has a flat power gain equal to X [dB].
%   Otherwise, ATYPE can be 'fixpower' if the amplifier takes the gain so
%   as to have an output average power for channel ceil(Nch/2) equal 
%   to X [mW], Nch being the number of channels.
%   This options works only with channels separated (see CREATE_FIELD).
%
%   AMPLIFLAT(X,ATYPE,OPTIONS) has the additional variable OPTIONS to
%   insert the amplified spontaneous emission (ASE) noise.
%   OPTIONS is a structure whose fields can be:
%
%   OPTIONS.f:     [dB] is the optical ASE noise figure, which corresponds 
%                  to a one-sided ASE power spectral density, on two 
%                  polarizations, N0 = F*(Gain-1)*h*nu, with Gain the 
%                  amplifier gain, h the Planck's constant and nu the 
%                  channel central frequency.
%                  Hence, ASE power on a frequency band B is Pase = N0*B.
%   OPTIONS.asepol:If it's 'asex' allows to force to zero the ASE noise 
%                  added to GSTATE.FIELDY, while for 'asey' allows to force
%                  to zero the noise added to GSTATE.FIELDX.
%   OPTIONS.noise: A matrix containing user's defined complex, unit
%                  variance and zero mean Gaussian ASE noise samples. 
%                  OPTIONS.noise must have the same size of 
%                  [GSTATE.FIELDX, GSTATE.FIELDY] and must be read in that 
%                  way.
%
%   If the amplifier does not generate ASE noise, don't set OPTIONS.
%   
%   Note: AMPLIFLAT assumes the same gain for both polarizations. 
%


global CONSTANTS;  % physical constants: a structure defined in reset_all.m
global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

HPLANCK = CONSTANTS.HPLANCK;      % Planck's constant [J*s]
CLIGHT = CONSTANTS.CLIGHT;      % speed of light [m/s]

[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 

atype = lower(atype);
switch atype
    case 'gain'
        gain = 10^(x*0.1);
    case 'fixpower'                 % gain -> out power = x for ch.1
        if nfc == GSTATE.NCH
            midch = ceil(nfc/2);
            avge = avg_power(midch,'abs'); % average energy per bit ch.midch
            gain = x/avge;    
        else
            error(['''fixpower'' works',...
                ' only for channels separated']);
        end
                                                    
    otherwise
        error('wrong string atype');
end
GSTATE.FIELDX = GSTATE.FIELDX*sqrt(gain);
isy = ~isempty(GSTATE.FIELDY);
if isy
    GSTATE.FIELDY = GSTATE.FIELDY*sqrt(gain);
end




        
       
