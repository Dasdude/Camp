function S1=wallop(w1,sdata,plotflag)
%WALLOP  Calculates (and plots) a Wallop spectral density.
% 
% CALL:  S = wallop(w,data,plotflag); 
%        S = wallop(wc,data,plotflag);
%
%        S    = a struct containing the spectral density, see datastructures.
%        w    = angular frequency (default linspace(0,3,257))
%        wc   = angular cutoff frequency (default 33/Tp)
%        data = [Hm0 Tp M]
%               Hm0 = significant wave height (default 7 (m))
%               Tp  = peak period (default 11 (sec))
%               M   = shape factor, i.e. slope for the high frequency
%                     part (default depending on Hm0 and Tp, see below)
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the spectrum.
%
%  The WALLOP spectrum parameterization used is 
%
%     S(w)=Bw*Hm0^2/wp*(wp./w).^M.*exp(-M/4*(wp./w).^4);
% where
%     Bw = normalization factor
%     M  = abs((log(2*pi^2)+2*log(Hm0/4)-2*log(Lp))/log(2));
%     Lp = wave length corresponding to the peak frequency, wp.
%
%  If M=5 it becomes the same as the JONSWAP spectrum with 
%  peak enhancement factor gamma=1 or the Pierson-Moskowitz spectrum. 
%
% Example: 
%   S = wallop(1.1,[6.5 10]), wspecplot(S)
%  
% See also  jonswap, torsethaugen, simpson

% References:
% Huang, N.E., Long, S.R., Tung, C.C, Yuen, Y. and Bilven, L.F. (1981)
% "A unified two parameter wave spectral model for a generous sea state"
% J. Fluid Mechanics, Vol.112, pp 203-224

% Tested on: matlab 6.0, 5.3
% History:
% revised jr 03.04.2001
%  - added wc to input 
%  - updated information
% revised pab 18.02.2000
%  - normalization so that int S(w) dw = m0
% revised pab 24.01.2000
%  - updated note to 'Wallop Hm0='....
% by pab 01.12.99

monitor=0;

if nargin<3|isempty(plotflag)
  plotflag=0;
end

% Old call  
%if nargin<1|isempty(w1)
%  w=linspace(0,3,257).';
%end

M=[];
if nargin<2|isempty(sdata)
  sdata=[7 11];
else
  switch length(sdata)
    case 1, sdata=[sdata 11];
    case 3, M=sdata(3); sdata=sdata(1:2);
  end
end %

if nargin<1|isempty(w1), wc = 33/sdata(2)
elseif length(w1)==1,    wc = w1; 
else w = w1 ; end
nw = 257;
if isempty(w), w = linspace(0,wc,nw).'; end
 
n=length(w);

S1=createspec;
S1.S=zeros(n,1);
S1.w=w;
S1.norm=0; % The spectrum is not normalized

Hm0=sdata(1);
Tp=sdata(2);
S1.note=['Wallop, Hm0 = ' num2str(Hm0)  ', Tp = ' num2str(Tp)];
wp=2*pi/Tp;

if monitor
  disp(['Hm0, Tp      = ' num2str([Hm0 Tp])])
end

if isempty(M),
  kp=w2k(wp,0,inf); % wavenumber at peak frequency
  Lp=2*pi/kp; % wave length at the peak frequency
  M=abs((log(2*pi^2)+2*log(Hm0/4)-2*log(Lp))/log(2));
end

%Bw=0.06238*M^((M-1)/4)/(4^((M-5)/4)*gamma((M-1)/4))*(1+0.7458*(M+2)^(-1.057));

Bw = M^((M-1)/4)/(4^((M-5)/4)*gamma((M-1)/4))/16;

% for w>0 % avoid division by zero
k=find(w>0);
S1.S(k)=Bw*Hm0^2/wp*(wp./w(k)).^M.*exp(-M/4*(wp./w(k)).^4);

if plotflag
  wspecplot(S1,plotflag)
end

