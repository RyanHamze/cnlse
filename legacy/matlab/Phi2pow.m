function pwr=phi2pow(phi,L,alpha,gam,G,nspan)

%PHI2POW Convert nonlinear phase into power.
%   PWR=PHI2POW(PHI,L,ALPHA,GAM,G,NSPAN) yields the transmitted power 
%   PWR [mW] corresponding to a cumulated nonlinear phase into the overall 
%   optical link equal to PHI [rad]. PHI is the self-phase modulation (SPM)
%   rotation induced by the link over a constant signal of power PWR.
%
%   L is a vector of size equal to the number of nonlinear fibers in the 
%   link and contains the fiber lengths in [m]. ALPHA are the fibers 
%   attenuation [dB/km], GAM are the fibers nonlinear coefficients, usually 
%   called gamma [1/W/m].
%   G is the net gain [dB] from fiber to fiber (see the example).
%
%   If L,ALPHA,GAM,G are of length 1, the same value is used for a total 
%   of NSPAN spans.
%   
%
%   E.g. Two-fiber periodic link of N spans. One of the NSPAN periods is:
%
%   ampli0     fiber 1     ampli1        fiber 2     ampli2=ampli0
%               ___                       ___    
%   |\         /   \       |\            /   \       |\
%   | \       |     |      | \          |     |      | \
%   |  \       \   /       |  \          \   /       |  \ 
%   |   >------------------|   >---------------------|   >--- ...
%   |  / 0               L1|  /                 L1+L2|  /
%   | /                    | /                       | /
%   |/                     |/                        |/ 
%   gain  attenuation a1   gain    attenuation a2    gain   
%   G0    nl-index gam1    G1      nl-index gam2     G2=G0
%
%
%   Call the function with:
%   L = [L1,L2], ALPHA = [a1,a2], GAM = [gam1,gam2], 
%   G=[G0,G1]
%
%   Note 1: If the laser is directly connected to the link, G0=0 dB
%   Note 2: transparent link of N equal spans -> L = L1, ALPHA=a1, 
%         G=0 and NSPAN = N.
%
%   The function returns PWR [mW] that entering the amplifier ampli0 gives
%   the cumulated nonlinear phase PHI.
%



%global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

Ll = length(L);
La = length(alpha);
Lgam = length(gam);
LG = length(G);

if (abs(Ll-La)+abs(Ll-Lgam)+abs(Lgam-LG)) ~= 0
    error('All fiber parameters must be of the same length');
end

alphalin = (log(10)*1e-4)*alpha;             % [m^-1]
indnz = find(alphalin);             % non-zero indexes
indz = find(alphalin == 0);         % zero indexes
Leff(indnz) = (1-exp(-alphalin(indnz).*L(indnz)))./alphalin(indnz);
Leff(indz) = L(indz);   % effective length [m]
loss = -alpha.*L*1e-3;       % fiber loss [dB]
netgain = [0,loss(1:La-1)]+G; % netgain(k)= net gain [dB] of (k-1)-th fiber + k-th ampli
cumgain = 10.^(cumsum(netgain)*0.1);    % cumulative gain: cumgain(k)*pwragv is the power input to the k-th fiber

pwr = phi/sum(Leff.*gam.*cumgain*nspan)*1e3;  % [mW]


