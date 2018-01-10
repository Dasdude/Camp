function Snew = specinterp(S,dt)
%SPECINTERP Interpolation and zero-padding of spectrum
%            to change Nyquist freq.
%
% CALL:  Snew = specinterp(S,dt)
%
%    Snew, S = spectrum structs (any type)
%         dt = wanted sampling interval (default as given by S, see spec2dt)
%              unit: [s] if frequency-spectrum, [m] if wave number spectrum  
%              
% To be used before simulation (e.g. spec2sdat) or evaluation of covariance
% function (spec2cov) to directly get wanted sampling interval.
% The input spectrum is interpolated and padded with zeros to reach
% the right max-frequency, w(end)=pi/dt, f(end)=1/2/dt, or k(end)=pi/dt.
% The objective is that output frequency grid should be at least as dense as
% the input grid, but the number of points in the spectrum is maximized to
% 2^13+1. 
% NB! Also zero-padding down to zero freq, if S does not start there. 
%     If empty input dt, this is the only effect.
%  
% See also  spec2cov, spec2sdat, covinterp, spec2dt

% Tested on:
% History:
% revised pab 12.10.2001
% -fixed a bug created 11.10.2001: ftype='k' now works OK
% revised pab 11.10.2001
% - added call to freqtype.m
% - fixed a bug: ftype=='f' now works correctly
% revised es 080600, revision of last revision, output matrix now OK
%  revised by es 22.05.00, output vectors columns, not rows
%  by es 13.01.2000, original ideas by sylvie, ir
  
Snew  = S;  
ftype = freqtype(S);
w     = getfield(S,ftype);
n     = length(w);

if strcmp(ftype,'f') %ftype==f
  dT=1/(2*w(n));% sampling interval=1/Fs
else % ftype == w og ftype == k
  dT=pi/w(n);
end
if nargin<2|isempty(dt),  dt=dT; end

% Find how many points that is needed
nfft   = 2^nextpow2(n-1);
dttest = dT*(n-1)/nfft;
while (dttest>dt) & (nfft<2^13)
 nfft=nfft*2;
 dttest=dT*(n-1)/nfft;
end;
nfft=nfft+1;

if strcmp(ftype,'w')
  Snew.w=linspace(0,pi/dt,nfft)';
elseif strcmp(ftype,'f') %freqtype==f
  Snew.f=linspace(0,1/2/dt,nfft)';
else
  Snew.k=linspace(0,pi/dt,nfft)';
end

if size(S.S,2)<2
  S.S=S.S'; % if vector, make it a row to match matrix (np x nf)
end
w=w(:);
dw=w(end)-w(end-1);

fmax = max(getfield(Snew,ftype));
if fmax>w(end)
  % add a zero just above old max-freq, and a zero at new max-freq
  % to get correct interpolation there
  w=[w;w(end)+dw; fmax ];
  S.S=[S.S zeros(size(S.S,1),2)];
end
if w(1)>0
  % add a zero at freq 0, and, if there is space, a zero just below min-freq
  if w(1)>dw
    w=[w(1)-dw;w];
    S.S=[zeros(size(S.S,1),1) S.S];
  end
  size(w)
  w=[0;w];
  S.S=[zeros(size(S.S,1),1) S.S];
end  

Snew.S=interp1q(w,S.S',getfield(Snew,ftype))';

if size(Snew.S,1)<2
  Snew.S=Snew.S.'; % if vector, make it a column
end




