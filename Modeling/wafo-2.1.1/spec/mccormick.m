function S1=mccormick(w1,sdata,plotflag)
%MCCORMICK Calculates (and plots) a McCormick spectral density.
% 
% CALL:  S = mccormick(w,data,plotflag); 
%        S = mccormick(wc,data,plotflag);
%
%        S    = a struct containing the spectral density, see datastructures.
%        w    = angular frequency (default linspace(0,33/Tp,257))
%        wc   = angular cutoff frequency (default 33/Tp)
%        data = [Hm0 Tp Tz]
%               Hm0  = significant wave height (default 7 (m))
%               Tp   = peak period (default 11 (sec))
%               Tz   = zero-down crossing period (default  0.8143*Tp)
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the spectrum.
%
%  The McCormick spectrum parameterization used is 
%
%     S(w) = (M+1)*(Hm0/4)^2/wp*(wp./w)^(M+1)*exp(-(M+1)/M*(wp/w)^M);
%  where
%     Tp/Tz=(1+1/M)^(1/M)/gamma(1+1/M)
%
% 
% Example: 
%   S = mccormick(1.1,[6.5 10]), wspecplot(S)
%
% See also  jonswap, torsethaugen, simpson


% References:
% M.E. McCormick (1999)
% "Application of the Generic Spectral Formula to Fetch-Limited Seas"
% Marine Technology Society, Vol 33, No. 3, pp 27-32


% Tested on: matlab 6.0, 5.3
% History:
% revised pab nov 2004
%  -replaced fmin with fminbnd   
% revised jr 03.04.2001
% - added wc to input 
% - updated information
% by pab 01.12.99


monitor=0;

if nargin<3|isempty(plotflag)
  plotflag=0;
end
if nargin<2|isempty(sdata)
  sdata=[7 11 11/1.228];
else
  switch length(sdata)
    case 1, sdata=[sdata, 11, 11/1.228];
    case 2, sdata=[sdata, 11/1.228];
  end
end 

w = [];
if nargin<1|isempty(w1), wc = 33/sdata(2);
elseif length(w1)==1,    wc = w1; 
else w = w1 ; end
nw = 257;
if isempty(w), w = linspace(0,wc,nw).'; end

% Old definition of w 
%if nargin<1|isempty(w)
%  w=linspace(0,33/sdata(2),257).';
%end


n=length(w);
S1=createspec;
S1.S=zeros(n,1);
S1.w=w;
S1.norm=0; % The spectrum is not normalized


Hm0=sdata(1);
if 0,
  Tz=sdata(2);
  Tp=sdata(3);
  S1.note=['McCormick, Hm0, Hm0 = ' num2str(Hm0)  ', Tm01 = ' num2str(Tp)];
else
  Tp=sdata(2);
  Tz=sdata(3);
  S1.note=['McCormick, Hm0 = ' num2str(Hm0)  ', Tp = ' num2str(Tp) ', Tz = ' num2str(Tz) ];
end
wp=2*pi/Tp;

if monitor
  disp(['Hm0, Tp, Tz     = ' num2str([Hm0 Tp Tz])])
end

mvrs=version;ix=find(mvrs=='.');
if str2num(mvrs(1:ix(2)-1))>5.2,
  M=fminbnd(['((1+1./x).^(1./x)/gamma(1+1./x)-' num2str(Tp/Tz) ').^2' ],1,7);
else
  M=fmin(['((1+1./x).^(1./x)/gamma(1+1./x)-' num2str(Tp/Tz) ').^2' ],1,7);
end



% for w>0 % avoid division by zero
k=find(w>0);
S1.S(k)=(M+1)*(Hm0/4)^2/wp*(wp./w(k)).^(M+1).*exp(-(M+1)/M*(wp./w(k)).^M);
             

if plotflag
  wspecplot(S1,plotflag)
end
