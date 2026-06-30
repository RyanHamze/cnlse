function zbrf=fiber(x,flag)

%FIBER Optical fiber in the nonlinear regime.
%   FIBER(X,FLAG) solves the nonlinear Schroedinger equation (NLSE) in 
%   absence of polarization effects, or the Coupled-NLSE (CNLSE) with 
%   polarization effects. X is a structure of fields:
%
%       X.length:     fiber length [m]
%       X.alphadB:    fiber attenuation [dB/km]
%       X.aeff:       fiber effective area [um^2]
%       X.n2:         fiber nonlinear index [m^2/W]
%       X.lambda:     wavelength [nm] of X.disp 
%       X.disp:       fiber chromatic dispersion coefficient [ps/nm/km]
%                     @X.lambda
%       X.slope:      fiber slope, i.e. derivative of X.disp [ps/nm^2/km]
%       X.dzmax:      max. step for the split-step algorithm [m]
%       X.dphimax:    max. nonlinear phase rotation in each step [rad]
%
%   For the solution of the CNLSE, i.e. with two polarizations, there are
%   also the following additional parameters:
%
%       X.dgd:        fiber average differential group delay [symbols]
%       X.nplates:    number of waveplates or trunks for PMD emulation
%       X.manakov:    'yes': Solve the Manakov equation. 'no': Solve the
%                     CNLSE. Default: 'no'.
%
%   In the general case with two polarizations the fiber is the 
%   concatenation of randomly oriented polarization mantaining fibers 
%   (PMF). The user can force the use of a single PMF by adding the 
%   following optional parameters:
%
%       X.db0:        birefringence of the PMF fiber at GSTATE.FN=0 
%       X.theta:      azimuth [rad] of the PMF fiber
%       X.epsilon:    ellipticity [rad] of the PMF fiber
%
%   The NLSE is solved by a split-step Fourier algorithm with a variable
%   step so as to have a maximum nonlinear phase rotation into each step
%   equal to X.dphimax. However, the step cannot be larger than X.dzmax. 
%
%   The CNLSE uses the same rules except that each X.length/X.nplates
%   meters a new PMF is generated.
%
%   Alternatively, the step can be chosen adaptively basing the choice on a
%   target local truncation error (NLSE only). In such a case the following 
%   parameters should be added to X:
%
%       X.ltol:       local truncation error, i.e. max distance between the
%                     field obtained by moving once or twice in a step.
%       X.dphiadapt:  true/false. True: the local truncation error method
%                     is applied only in the first step and used to correct 
%                     X.dphimax. After the first step the SSFM proceeds  
%                     using the approach based on X.dphimax. Default: false
%
%   FLAG is a string of four characters governing the type of propagation. 
%   The first character is 'g' if GVD (i.e. beta2,beta3) is on or '-' in 
%   absence of GVD. Note that with 'sepfields' in reset_all.m this function 
%   accounts for the walkoff effect even with the GVD flag set to '-'. 
%   The second character is 'p' for propagation of a polarized field in 
%   presence of birefringence and PMD or '-' in absence of such effects. 
%   The third is 's' if SPM is on or '-' in absence of SPM. 
%   Likewise, the fourth character is 'x' or '-' in presence/absence of 
%   XPM. The most complete case is FLAG='gpsx' and corresponds to 
%   propagation in presence of fiber GVD+PMD+SPM+XPM. 
%
%   E.g. Propagation with GVD+XPM, without PMD and SPM -> FLAG='g--x'
%
%   The fourth character of FLAG is active only with channels separated 
%   (see option 'sepfields' in CREATE_FIELD). In this case, the propagation 
%   neglects the effect of four-wave mixing, which can be taken in account 
%   only by combining all channels into a unique field and hence it is a 
%   special case of SPM.
%
%   OUT=FIBER(X,FLAG) returns in OUT a struct containing the birefringence
%   parameters used by fiber:
%   
%   OUT.db0 = birefringence [rad] at GSTATE.FN=0 (see RESET_ALL).
%   OUT.theta = azimuth [rad] of all the PMFs composing the fiber.
%   OUT.epsilon = ellepticity [rad] of all the PMFs composing the fiber.
%   OUT.dgd = DGD [symbols].
%   OUT.lcorr = length [m] of each PMF trunk.
%   OUT.betat = beta(omega), i.e. scalar phase shift [rad] including GVD,
%       slope,etc, where omega/2/pi is the vector of FFT frequencies. betat
%       is common to both polarizations.
%   OUT.db1 = differential phase shift [rad] induced by PMD.
%
%   OUT can be used to recover the PMD transfer matrix of the fiber (see 
%   INVERSE_PMD).
%
%   The CNLSE is described as the concatenation of X.nplates PMF trunks,
%   each with principal states of polarization randomly distributed over the 
%   Poincare sphere, if not otherwise specified. Each PMF has constant DGD 
%   and randomly distributed birefringence. The nonlinearity is inserted 
%   after a certain number of trunks, depending on X.dzmax and X.dphimax. 
%   The diagram is the following:
%
%    ------------ -----------     --------- ------------ -----------
%   |  Lin step  | Lin step  |   | NL STEP |  Lin step  | Lin step  |
%   |    PMF 1   |   PMF 2   |...|         |    PMF k   |   PMF k+1 |...
%    ------------ -----------     --------- ------------ -----------
%
%   where PMF k, k=1,2,... is a PMF fiber randomly chosen on the Poincare
%   sphere.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global CONSTANTS;  % CONSTANTS is a global structure variable.
CLIGHT = CONSTANTS.CLIGHT;      % speed of lightin vacuum [m/s]
global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

SAFETYFCT = 0.9;    % Safety step-reduction factor
DEF_PLATES = 100;   % number of waveplates when birefringence is on
[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
Nfft = GSTATE.NSYMB*GSTATE.NT;
nfnames = fieldnames(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK AND INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    error('Missing propagation type');
end
if ~any(strcmp(nfnames,'dzmax')) || (x.dzmax > x.length)
    x.dzmax = x.length; 
end

if any(strcmp(nfnames,'ltol')), % adaptive step-size
    if ~any(strcmp(nfnames,'dphimax')), x.dphimax = Inf;end
    trg.err = x.ltol;
    trg.safety = SAFETYFCT;
    if any(strcmp(nfnames,'dphiadapt')) && x.dphiadapt
        tolflag = 1; % adaptive first step, then constant phase method
    else
        tolflag = 2; % adaptive steps with target local error
    end
else
    trg = [];
    tolflag = 0; % constant phase method
end

fls = [0 0 0 0];  % flag converted into array of numbers
switch lower(flag)
    case '----'  % only attenuation
        dphimaxt = Inf;
        dzmaxt = x.length;
    case 'g---'  % only GVD
        fls(1) = 1;
        dphimaxt = Inf;
        dzmaxt = x.length;
    case '-p--'  % only PMD
        fls(2) = 1;
        dphimaxt = Inf;
        dzmaxt =x.length;        
    case '--s-'  % only SPM
        fls(3) = 1;
        if nfc == 1
            dphimaxt = Inf; % with one field we have the exact solution
            dzmaxt = x.length;
        else
            dphimaxt = x.dphimax;
            dzmaxt = x.dzmax;     
        end
    case '---x'  % only XPM 
        if nfc == 1
            error('flag ''---x'' available only for channels separated');
        end
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;              
        fls(4) = 1;
    case 'gp--'  % GVD+PMD        
        fls(1:2) = 1;
        dphimaxt = Inf;
        dzmaxt = x.length;                       
    case 'g-s-'  % GVD+SPM
        fls([1 3]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;              
    case 'g--x'  % GVD+XPM
        if nfc == 1
            error('flag ''g--x'' available only for channels separated');
        end        
        fls([1 4]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;    
    case '-ps-'  % PMD+SPM
        fls([2 3]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;              
    case '-p-x'  % PMD+XPM
        if nfc == 1
            error('flag ''-p-x'' available only for channels separated');
        end        
        fls([2 4]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;              
    case '--sx'  % SPM+XPM
        fls(3) = 1;
        if (nfc ~= 1)
            fls(4) = 1;
            dphimaxt = x.dphimax;
            dzmaxt = x.dzmax; 
        else
            dphimaxt = Inf;
            dzmaxt = x.length;
        end
    case 'g-sx'  % GVD+SPM+XPM
        fls([1 3]) = 1;
        if (nfc ~= 1), fls(4) = 1;end;        
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;     
    case '-psx'  % PMD+SPM+XPM
        fls(2:3) = 1;
        if (nfc ~= 1), fls(4) = 1;end;        
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;
    case 'gps-'  % GVD+PMD+SPM
        fls(1:3) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;     
    case 'gp-x'  % GVD+PMD+XPM
        if nfc == 1
            error('flag ''gp-x'' available only for channels separated');
        end     
        fls([1 2 4]) = 1;
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;             
    case 'gpsx'  % GVD+PMD+SPM+XPM
        fls(1:3) = 1;
        if (nfc ~= 1), fls(4) = 1;end;        
        dphimaxt = x.dphimax;
        dzmaxt = x.dzmax;     
        
    otherwise
        error('wrong flag. E.g. flag can be ''g---'',''gp--'',''-s--'', etc');
end

isy = ~isempty(GSTATE.FIELDY); 
isv = (fls(2) == 1) || isy;
if fls(2) == 1
    if ~any(strcmp(nfnames,'manakov'))
        x.manakov = 'no';
    end 
    if ~any(strcmp(nfnames,'dgd')), error('Missing DGD in fiber');end
    isdb0 = any(strcmp(nfnames,'db0'));         % exist db0?
    istheta = any(strcmp(nfnames,'theta'));
    isepsilon = any(strcmp(nfnames,'epsilon'));
    ispmf = isdb0+istheta+isepsilon;
    if ispmf == 3   % PMF fiber
        x.nplates = length(x.theta);
        brf.db0 = x.db0;          % delta beta 0 of a PMF fiber
        brf.theta=x.theta;        % theta and epsilon are the angles of the PMF
        brf.epsilon=x.epsilon;
        dgdrms = x.dgd/x.nplates; % DGD per trunk
    elseif ispmf == 0 
        if ~any(strcmp(nfnames,'nplates'))
            x.nplates = DEF_PLATES;    % default value
        end
        brf.db0 = rand(x.nplates,1)*2*pi - pi;   % general case: No PMF.   
        brf.theta=rand(x.nplates,1)*pi - 0.5*pi;     % azimuth: uniform R.V.;   
        brf.epsilon=0.5*asin(rand(x.nplates,1)*2-1); % uniform R.V. over the Poincare sphere       
        dgdrms = sqrt((3*pi)/8)*x.dgd/sqrt(x.nplates);
        % DGD rms per trunk. x.dgd is the average value. 3pi/8 comes from
        % the Maxwellian distribution of DGD
    else                                             
        error('Missing one of db0, theta or epsilon in fiber');
    end
    brf.dgd = x.dgd;
    dgdrms = dgdrms/GSTATE.SYMBOLRATE; % convert in [ns]
    if ~isy, 
        GSTATE.FIELDY = zeros(nfr,nfc); % create y component
        warning('optilux:fiber',['You are working with two polarizations ',...
            'but in create_field you initialized just one polarization']);
    end
else
    dgdrms = 0; % turn off all polarization and birefringence effects
    brf.db0 = 0;
    brf.theta = 0;
    brf.epsilon = 0;
    ispmf = 1;
    x.manakov = 'no';
    x.nplates = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alphalin = (log(10)*1e-4)*x.alphadB;      % [m^-1]
if alphalin == 0
    Leff = x.length;      % effective length [m]
else
    Leff = (1-exp(-alphalin*x.length))/alphalin;
end
b20 = -x.lambda^2/2/pi/CLIGHT*x.disp*1e-6; % beta2 [ns^2/m] @ lambda
b30 = (x.lambda/2/pi/CLIGHT)^2*(2*x.lambda*x.disp+x.lambda^2*x.slope)*1e-6;  
                                     % beta3 [ns^3/m] @ lambda
b30 = b30*fls(1);   % insert b30 only with GVD flag on.                                     

maxl = max(GSTATE.LAMBDA);
minl = min(GSTATE.LAMBDA); 
lamc = 2*maxl*minl/(maxl+minl); %central wavelength: 1/lamc = 0.5(1/maxl+1/minl)

% Domega_ik: [1/ns]. "i" -> at ch. i, "0" -> at lambda, "c" -> at lamc
Domega_i0 = 2*pi*CLIGHT*(1./GSTATE.LAMBDA-1/x.lambda);    
Domega_ic = 2*pi*CLIGHT*(1./GSTATE.LAMBDA-1/lamc);    
Domega_c0 = 2*pi*CLIGHT*(1./lamc-1/x.lambda); 
b1 = b20*Domega_ic+0.5*b30*(Domega_i0.^2-Domega_c0.^2);  %ch's beta1 [ns/m]
if nfc == 1
    beta1 = 0;    % [ns/m] @ lamc
    Domega_i0 = 2*pi*CLIGHT*(1./lamc-1/x.lambda);        
	gam = 2*pi*x.n2./(lamc*x.aeff)*1e18; % nonlinear index [1/mW/m]
else                                     
    beta1 = b1;   % [ns/m] @ GSTATE.LAMBDA
	gam = 2*pi*x.n2./(GSTATE.LAMBDA*x.aeff)*1e18; % nonlinear index [1/mW/m]
end
beta2 = b20+b30*Domega_i0;  % beta2 [ns^2/m]@ lamc (nfc=1) or GSTATE.LAMBDA
                            % (nfc ~= 1).
beta2 = beta2*fls(1);    % insert GVD only if the GVD flag is on. 
% Note:  fls(1)=0 means that we are neglecting local GVD, but we account
% for the walk-off effect.

Dch = x.disp+x.slope*(GSTATE.LAMBDA-x.lambda);  % dispersion of the channels
% Note: channel 1 is the reference channel. All the other channels shift
%   (beta1 effect) with respect to it.

Ld = Inf*ones(1,GSTATE.NCH);
fzb2 = find(Dch);
Ld(fzb2) = 1./(GSTATE.SYMBOLRATE^2*abs(x.lambda^2/2/pi/CLIGHT*Dch(fzb2)*1e-6));
% dispersion length [m]
if b30 ~= 0
    Lds = 1./(GSTATE.SYMBOLRATE^3*abs(b30));      % slope length [m]
else
    Lds = Inf;
end

betat = zeros(Nfft,nfc);    % beta, whose Taylor coeffs are the betai
db1 = zeros(Nfft,nfc);  % delta beta, i.e. phase shift of DGD.
omega = 2*pi*GSTATE.SYMBOLRATE*GSTATE.FN';     % angular frequency [rad/ns]

for kch=1:nfc
    betat(:,kch) = omega*beta1(kch)+0.5*omega.^2*beta2(kch)+...
        omega.^3*b30/6;   % beta coefficient [1/m]
    if fls(2) == 1
        db1(:,kch) = dgdrms*omega;  % DGD phase shift  
        % since the step is non-uniform, the DGD will be normalized in
        % matrix_step
    end
end

Lnl = 1./(gam.*GSTATE.POWER);    % nonlinear length [m]

%%%%%%%%%%%%%%%%%% UPDATE DELAY AND CUMULATED DISPERSION %%%%%%%%%%%%%%%%%%
loc_delay = x.length*GSTATE.SYMBOLRATE.*b1;
GSTATE.DELAY = GSTATE.DELAY+ones(1+isy,1)*loc_delay;
GSTATE.DISP = GSTATE.DISP + ones(1+isy,1)*fls(1)*Dch*x.length*1e-3;% cum .dispersion [ps/nm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SSFM PROPAGATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%

if tolflag == 2 % adaptive step, local truncation error method
    if isv  % with polarization
        error('adaptive step available in absence of polarization effects');
    else    % scalar propagation
        [firstdz,ncycle,GSTATE.FIELDX]=scalar_a_ssfm(GSTATE.FIELDX,betat,...
            dzmaxt,dphimaxt,gam,alphalin,Nfft,nfc,x.length,trg,fls);
    end
else            % constant phase x step
    if isv  % with polarization
        [firstdz,ncycle,GSTATE.FIELDX,GSTATE.FIELDY,brf]=matrix_ssfm(...
            GSTATE.FIELDX,GSTATE.FIELDY,betat,db1,dzmaxt,dphimaxt,gam,...
            alphalin,nfc,x.length,x.nplates,x.manakov,fls,brf);
        if nargout, zbrf=brf;end
    else    % scalar propagation
        [firstdz,ncycle,GSTATE.FIELDX]=scalar_ssfm(GSTATE.FIELDX,betat,...
            dzmaxt,dphimaxt,gam,alphalin,Nfft,nfc,x.length,fls,tolflag,trg);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRINT SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
function [firstdz,ncycle,u]=scalar_ssfm(u,betat,dzmaxt,dphimaxt,gam,...
    alphalin,Nfft,nfc,Lf,fls,tolflag,trg)

%SCALAR_SSFM SSFM algorithm for the scalar NLSE 
%   [FIRSTDZ,NCYCLE,U]=SCALAR_SSFM(U,BETAT,DZMAXT,DPHIMAXT,GAM,ALPHALIN,...
%   NFFT,NFC,LF,FLS,TOLFLAG,TRG)  
%   solves the Nonlinear Schrodinger equation through the split step Fourier 
%   method (SSFM). This function works for a scalar field U.
%
%   BETAT is the the (scalar) beta(omega) coefficent. DZMAXT is the
%   largest step allowed by nonlinear effects. DPHIMAXT is the largest
%   nonlinear phase rotation x step. GAM is the nonlinear coefficient.
%   ALPHALIN is the attenuation [m^-1]. NFFT is the number of FFT points.
%   NFC is the number of fields. LF is the fiber length [m]. 
%   FLS is the propagation flag. TOLFLAG=1 corrects DPHIMAXT adaptively.
%   TRG is a struct containing the flags of the adaptive step algorithm.
%
%   This function operates over the global field GSTATE.FIELDX in absence 
%   of PMD effects. In presence of PMD see MATRIX_SSFM.
%
%   The function uses the symmetric SSFM if the step is chosen adaptively 
%   on the basis of the local truncation error. Otherwise a basic SSFM is
%   used.
%
%   The function returns the electric field U after propagation, the
%   first step FIRSTDZ and the number of cycles NCYCLE used in the SSFM.

gamrep = repmat(gam,Nfft,1);
dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,u,[],false);
halfalpha = 0.5*alphalin; 
ncycle = 1;      % number of cycles
if tolflag == 1 % adaptive search of dphimax
    if dz >= dzmaxt
        Umax = max(real(u).^2+imag(u).^2);
        maxpow = max(gam.*Umax);
        if alphalin == 0
            dphimaxt = maxpow*dz;
        else
            dphimaxt = maxpow*(1-exp(-alphalin*dz))/alphalin;
        end
    end
    dzini = dz;
    zdone=0;
    while zdone == 0
        [u ncycle nrej zdone dz]=adaptssfm(u,zdone,dz,alphalin,...
            gamrep,nfc,fls,betat,halfalpha,trg,0,0);
    end
    if dz > dzmaxt, dz = dzmaxt;end
    % on the basis of the first step correct the max phase rotation x step
%     fprintf('dphiINI=%.5f\n',dphimaxt)
    dphimaxt = dphimaxt*(1-exp(-alphalin*zdone))/(1-exp(-alphalin*dzini));
%      fprintf('dphiADAPT= %.5f\n',dphimaxt); 
    firstdz = zdone;
    zprop = zdone+dz;
    ncycle = ncycle+1;
else
    firstdz = dz;
    zprop = dz;             % running distance [m]
end
while zprop < Lf        % all steps except the last
    u = nl_step(alphalin,gamrep,dz,u,nfc,fls(3),fls(4));
                                                     % 1/3) NON-linear step
    
    u = lin_step(betat*dz,u); % 2/3) linear step    
    
    u = u*exp(-halfalpha*dz);    % 3/3) attenuation
    
    dz = nextstep(dzmaxt,dphimaxt,gam,alphalin,u,[],false);
    zprop = zprop+dz;
%     semilogy(ncycle,dz,'b+')
    ncycle = ncycle+1;
%     fprintf('zdone = %.3f\t  dz=%.3f\t \t  ncycle=%d\n',zprop,dz,ncycle)

end 
last_step = Lf-zprop+dz;
u = nl_step(alphalin,gamrep,last_step,u,nfc,fls(3),fls(4)); 
                                                          % NON-linear step
u = lin_step(betat*last_step,u); % linear step    

u = u*exp(-halfalpha*last_step);    % attenuation


%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function dz_nl=nextstep(dzmax,phimax,gam,alphalin,ux,uy,isv)

%NEXTSTEP step for the SSFM algorithm
%   DZ=NEXTSTEP(DZMAX,PHIMAX,GAM,ALPHALIN,U,NFC,ISV) evaluates the step of  
%   the SSFM algorithm for a fiber having nonlinear gamma coefficient GAM 
%   [1/mW/km] and attenuation ALPHALIN [m^-1]. U is the matrix of electrical
%   fields. NFC is the number of channels (equal to 1 if all channels are 
%   combined into a unique field). ISV=1 if the y component is on.
%   The step corresponds to a maximum nonlinear phase rotation equal to 
%   PHIMAX, under the constraint that it cannot be greater than DZMAX.

if isv % Umax is the largest normalized power
    Umax = max(real(ux).^2+imag(ux).^2+real(uy).^2+imag(uy).^2); 
else
    Umax = max(real(ux).^2+imag(ux).^2); 
end    
Pmax = max(gam.*Umax);   % largest gamma*power
leff = phimax/Pmax;             % effective length of the step
dl = alphalin*leff;             % ratio effective length/attenuation length

if dl >= 1
    dz_nl = dzmax;
else
    if alphalin == 0
        step = leff;
    else
        step = -1/alphalin*log(1-dl);
    end
    if step > dzmax
        dz_nl = dzmax;
    else
        dz_nl = step;
    end
end

%--------------------------------------------------------------------------



%-----------------------------------------------------------
function u = lin_step(betaxdz,u)

%LIN_STEP linear fiber in scalar propagation
%   U=LIN_STEP(BETAXDZ,U) propagates the electric field U into a purely
%   linear fiber having beta*dz factor BETAXDZ, i.e. a fiber
%   with transfer function exp(-i*beta(omega)*dz), where  betat is
%   approximated by its Taylor series in omega.


Hf = fastexp(-betaxdz); %  Fast exponential: fastexp(x) = exp(i*x)

u = ifft( fft(u) .* Hf);

%-----------------------------------------------------------
function u = nl_step(alphalin,gam,dz,u,nfc,spm,xpm)

%NL_STEP nonlinear fiber in scalar propagation
%   U=NL_STEP(ALPHALIN,GAM,DZ,U,NFC,SPM,XPM) propagates the electric field 
%   U into a purely nonlinear fiber having a nonlinear fiber having
%   nonlinear coefficient GAM [1/mW/km], attenuation ALPHA_LIN [m^-1]. NFC 
%   is the number of channels. SPM and XPM are 1 or 0 if self-phase 
%   modulation and cross-phase modulation have to be taken in account, 
%   respectively.


if alphalin == 0
    leff = dz;
else
    leff = (1-exp(-alphalin*dz))/alphalin;  % effective length of dz
end
pow = (real(u).^2+imag(u).^2);
if xpm
    if spm
        pow = 2*sum(pow,2)*ones(1,nfc)-pow; % sum(.,2): row sum
    else
        pow = 2*(sum(pow,2)*ones(1,nfc)-pow);        
    end
else
    if ~spm
        return
    end
end
u = u.*fastexp(-gam.*pow*leff);



%-----------------------------------------------------------

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
