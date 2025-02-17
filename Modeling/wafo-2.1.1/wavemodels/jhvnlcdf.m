function f = jhvnlcdf(Hd,Vcf,Hm0,Tp,gam,tail)
%JHVNLCDF Joint (Vcf,Hd) CDF for non-linear waves with JONSWAP spectrum.
%
%  CALL: f = jhvnlcdf(Hd,Vcf,Hm0,Tp,Gamma,tail)
% 
%   f   = CDF evaluated at (Vcf,Hd)
%   Hd  = zero down crossing wave height [m] 
%   Vcf = crest front velocity    [m/s]
%   Hm0 = significant wave height [m]
%   Tp  = Spectral peak period    [s]
% Gamma = Peakedness parameter of the JONSWAP spectrum  
%  tail = 1 if upper tail is calculated   
%         0 if lower tail is calulated (default)
%  
% JHVNLCDF approximates the joint CDF of (Vcf, Hd), i.e., crest front
% velocity (Ac/Tcf) and wave height, for a Gaussian process with a
% JONSWAP spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Vcf,Hd) data for 13 classes of
% GAMMA between 1 and 7. Between 50000 and 150000 zero-downcrossing waves were
% simulated for each class of GAMMA.
% JHVNLCDF is restricted to the following range for GAMMA: 
%  1 <= GAMMA <= 7 
%
% Example:
% Hm0 = 6;Tp = 8; gam = 3.5;
% vc = 3;
% hc = 3;
% lowerTail = 0;
% upperTail = ~lowerTail  
% jhvnlcdf(hc,vc,Hm0,Tp,gam)           % Prob(Hd<Hc,Vcf<Vc)
% jhvnlcdf(hc,vc,Hm0,Tp,gam,upperTail) % Prob(Hd>Hc,Vcf>Vc)  
%  
%  % Conditional probability of steep and high waves given seastates
%  % i.e., Prob(Hd>hc,Vcf>vc|Hs,Tp)  
%  upperTail = 1;
%  Hs = linspace(2.5,11.5,10);
%  Tp = linspace(4.5,19.5,16);
%  [T,H] = meshgrid(Tp,Hs); 
%  pl = jhvcdf(hc,vc,H,T,gam,upperTail); 
%  p = jhvnlcdf(hc,vc,H,T,gam,upperTail);
%  v = 10.^(-6:-1);  
%  contourf(Tp,Hs,log10(p),log10(v))
%  xlabel('Tp'), ylabel('Hs'),  fcolorbar(log10(v))  
%  figure(2)  
%  contourf(Tp,Hs,log10(pl),log10(v))
%  xlabel('Tp'), ylabel('Hs'),  fcolorbar(log10(v))  
%    
%  
% See also  thvpdf

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.    
  
  
% History
% revised pab 09.08.2003
% changed input  
% validated 20.11.2002
% By pab 20.12.2000


error(nargchk(2,6,nargin))  

if (nargin < 6|isempty(tail)),tail = 0; end
if (nargin < 4|isempty(Tp)),Tp = 8; end
if (nargin < 3|isempty(Hm0)), Hm0 = 6; end
if (nargin < 5|isempty(gam)),
   gam = getjonswappeakedness(Hm0,Tp);
end

multipleSeaStates = any(prod(size(Hm0))>1|...
			prod(size(Tp))>1|...
			prod(size(gam))>1);
if multipleSeaStates
  [errorcode, Vcf,Hd,Hm0,Tp,gam] = comnsize(Vcf,Hd,Hm0,Tp,gam);
else
  [errorcode, Vcf,Hd] = comnsize(Vcf,Hd);
end
if errorcode > 0
  error('Requires non-scalar arguments to match in size.');
end
displayWarning = 0;
if displayWarning
  if any(any(Tp>5*sqrt(Hm0) | Tp<3.6*sqrt(Hm0)))
    disp('Warning: Hm0,Tp is outside the JONSWAP range')
    disp('The validity of the parameters returned are questionable')
  end
end
%dev = 2e-5;
c1 =[ 0.16183666835624   1.53691936441548   1.55852759524555];
c2 =[ 0.15659478203944   1.15736959972513   1];
Tm02 = Tp.*(polyval(c2,gam)./polyval(c1,gam));

Hrms = Hm0/sqrt(2);
Vrms = 2*Hm0./Tm02; % Erms

v = Vcf./Vrms;

hMax = 5;
eps2 = 1e-6;
normalizedInput = 1;


h = min(Hd./Hrms,hMax);

f = zeros(size(Hd));
% Only compute
k0 = find((Tp<5*sqrt(Hm0)) & (3.6*sqrt(Hm0)<Tp));
if any(k0)
  if multipleSeaStates
    h = h(k0);
    v = v(k0);
    Hm0 = Hm0(k0);
    Tp = Tp(k0);
    gam = gam(k0);    
  else
    k0 = 1:prod(size(Hd));
  end
  
  
  utprb = gaussq('jhvnlpdf',hMax,2*hMax,eps2/2,[],mean(v(:)),mean(Hm0(:)),mean(Tp(:)),mean(gam(:)),normalizedInput,7);
  if eps2<utprb
    warning('Check the accuracy of integration!')
  end
  
  if 0
    % This is a trick to get the html documentation correct.
    k = jhvnlpdf(1,1,2,3);
  end
  hlim  = h;

  lowerTail = 0;
  if tail==lowerTail,
    k       = find(h>1.3*v);
    hlim(k) = 1.3*v(k);
    
    f(k0) = gaussq('jhvnlpdf',0,hlim,eps2/2,[],v,Hm0,Tp,gam,normalizedInput,5)...
	    + gaussq('jhvnlpdf',hlim,h,eps2/2,[],v,Hm0,Tp,gam,normalizedInput,5); 
  else % upper tail
    k       = find(h<1.3*v);
    hlim(k) = 1.3*v(k);
    f(k0) = gaussq('jhvnlpdf',h,hlim,eps2/2,[],v,Hm0,Tp,gam,normalizedInput,7)...
      + gaussq('jhvnlpdf',hlim,hMax,eps2/2,[],v,Hm0,Tp,gam,normalizedInput,7); 
  end
end
return

