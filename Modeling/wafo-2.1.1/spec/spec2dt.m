function dT = spec2dt(S);
% SPEC2DT  Computes sampling interval from Nyquist frequency in spectrum
%
%  CALL : dT = spec2dt(S)
%
%         dT = sampling interval, unit: [m] if wave number spectrum, [s] else  
%         S  = spectrum struct  
%
%  Let wm be maximum frequency/wave number in spectrum,
%  then dT=pi/wm if angular frequency, dT=1/(2*wm) if natural frequency (Hz)

% Tested on Matlab 5.3
% Revised by es 25.05.00: help text + call to freqtype
  
ftype=freqtype(S);

n=length(getfield(S,ftype));

if strcmp(ftype,'w')
  dT=pi/(S.w(n));
elseif strcmp(ftype,'f') 
  dT=1/(2*S.f(n)); % sampling interval=1/Fs
else
    dT=pi/(S.k(n));
end

