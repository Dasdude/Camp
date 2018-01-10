function Sf=spa2time(S,w,theta,g),
% SPA2TIME Transform of spectrum from wave no. to frequency (used in spec2spec)
%
% CALL:  Sf = spa2time(S,w,theta,g)
%
%   Sf = frequency or directional spectrum, which depends on input
%   S  = spectrum struct
%  w,theta = lag of output (default depending on input spectrum)
%   g  = constant of gravity (default see gravity)
%
% Transforms a wave number spectrum to a frequency or directional spectrum 
% Input 1D gives output 1D, input 2D gives output 2D
% Rotation, transformation etc is preserved
%
% See also definitions, spec2spec

% Tested on Matlab 5.3
% History: revised by es 29.11.1999: norm. in both x and y (length(g) 1 OR 2)
%          by es 99.08.17

if nargin<3
  if isfield(S,'g')
    g=S.g;
  else
    g=gravity;
  end
end

if ~(strcmpi(S.type(end-2:end),'k1d')|strcmpi(S.type(end-2:end),'k2d'))
  error('Spectrum already in time domain')
end
Sf=S;
k=S.k(:)';
Sk=S.S;
Sf=rmfield(Sf,'k');
%if isfield(Sf,'phi')
%  Sf=rmfield(Sf,'phi'); % a rotation has no sense in freq-spectrum
%end
if nargin<2|isempty(w)
  w=linspace(0,k2w(k(end),0,S.h,g),length(k))';
else
  if length(w)==1 % then interpret w as dw (step length)
    w=(0:length(k)-1)*w;
  end
end
w=w(:)'; % row vector
Sf.w=w;

if strcmpi(S.type(end-2:end),'k1d')

  kw=w2k(w,0,S.h,g(1));
  in=(kw>k(1))&(kw<=k(end));
  Gf=zeros(size(k));
  Gf(in)=interp1(k,Sk,kw(in));
  Sf.S=zeros(size(w));
  Sf.S=Gf.*dkdw(w,kw,S.h,g(1));
  Sf.type='freq';
else
  if S.h<inf
    disp('This transformation for finite depth is not available yet')
    return
  end
  if nargin<3|isempty(theta)
    if S.k(1)<0
      theta=linspace(-pi,pi,2^7+1);
    else
      theta=linspace(-pi/2,pi/2,2^6+1);
    end
  end
  theta=theta(:);
  Sf.theta=theta;
  W=ones(size(theta))*w;
  TH=theta*ones(size(w));
  KW=W.^2;
  K=KW.*cos(TH)/g(1);  % non-equidistant K derived from (W,TH)
  K2=KW.*sin(TH)/g(end); % non-equidistant K2 derived from (W,TH)
  Gk=zeros(size(W));
  in=(K>k(1))&(K<=k(end))&(K2>S.k2(1))&(K2<=S.k2(end));
  Gk(in)=interp2(k,S.k2,Sk,K(in),K2(in),'*linear');
  Sf.S=Gk.*W.^3*2/g(1)/g(end);  % directional spectrum
  Sf=rmfield(Sf,'k2');
  if strcmpi(S.type,'rotk2d')
    Sf.type='rotdir';
  else
    Sf.type='dir';
  end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dk=dkdw(w,k,h,g)

if h==inf
  dk=2*w/g;
else
  dk=zeros(size(w));
  dk(k>0)=2*w(k>0)./(g*tanh(k(k>0)*h)+g*h*k(k>0).*(1-tanh(k(k>0)*h).^2));
end
return
