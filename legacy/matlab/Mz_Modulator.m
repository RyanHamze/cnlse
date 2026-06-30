function Eout=mz_modulator( Ein, modsig, options)
%MZ_MODULATOR modulate the optical field with a Mach-Zehnder Interferometer
%   E=MZ_MODULATOR(E,MODSIG, OPTIONS) modulates the optical field E using a 
%   Mach-Zehnder interferometer [1]. 
%   The parameter MODSIG is the electrical driving signal produced by 
%   ELECTRICSOURCE. OPTIONS is an optional structure whose fields can be:
%
%       - exratio : extinction ratio [dB] (default = inf)
%       - bias  : bias of the modulator (default = 0)
%       - amplitude: Vpi of the modulator (default = 1)
%       - nochirp: reduce effect of chirp due to finite exratio [2]
%         (default = false)
%   
%   See also LASERSOURCE, ELECTRICSOURCE, 



bias      = 0;
amplitude = 1;
exratio     = inf;
nochirp   = false;

if exist('options','var')
    checkfields(options,{'bias','amplitude','exratio','nochirp'});
    if isfield(options,'bias');
        bias = options.bias;
    end
    if isfield(options,'amplitude');
        amplitude = options.amplitude;
    end
    if isfield(options,'exratio');
        exratio = options.exratio;
    end  
    if isfield(options,'nochirp');
        nochirp = options.nochirp;
    end      
end

exr_lin = 10^(-exratio/10);
gamma = (1-sqrt(exr_lin))/(sqrt(exr_lin)+1);

if nochirp
    Phi_U = (modsig *          pi/(1+gamma^2) * amplitude + bias);
    Phi_L = (modsig * -gamma^2*pi/(1+gamma^2) * amplitude + bias);
else
    Phi_U = (modsig *  pi/2 * amplitude + bias);
    Phi_L = (modsig * -pi/2 * amplitude + bias);
end

% Any phase shift due to the interferometric structure is neglected for the
% sake of simplicity
Eout  = j*Ein.*(fastexp(Phi_L)-gamma*fastexp(Phi_U))/(1+gamma);
