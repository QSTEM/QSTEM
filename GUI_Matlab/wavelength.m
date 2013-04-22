function lambda = wavelength(v0)
% returns the relativistic electron wavelength for an input voltage in keV

emass = 510.99906;   % electron rest mass in keV
hc = 12.3984244;     % h*c

lambda = hc./sqrt(v0.*(2*emass+v0));    % in Angstroem

