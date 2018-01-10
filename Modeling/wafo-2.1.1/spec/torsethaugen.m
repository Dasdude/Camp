function [S, Sw, Ss]=torsethaugen(w1,sdata,plotflag)
% TORSETHAUGEN Calculates a double peaked (swell + wind) spectrum 
%
%  CALL: [S, Ss, Sw] = torsethaugen(w,data,plotflag); 
%        [S, Ss, Sw] = torsethaugen(wc,data,plotflag); 
%
%   S, Ss, Sw = spectrum struct for total, swell and wind, respectively.  
%        w    = angular frequency (default linspace(0,wc,257))
%        wc   = angular cutoff frequency (default 33/Tp)
%        data = [Hm0 Tp A]
%               Hm0 = Significant wave height      (default 7)
%               Tp  = 2*pi/wp, primary peak period (deafult 11)
%               A   = alpha, normalization factor, (default -1)
%                   A<0  : A calculated by integration so that int S dw =Hm0^2/16
%                   A==0 : A = (1+f1(N,M)*log(gammai)^f2(N,M))/gammai, original 
%                         parametric normalization  
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the spectrum.
%
% For zero values, NaN's or parameters not specified in DATA the
% default values are used. 
%
% Model:  S(w)=Ss(w)+Sw(w) 
% where Ss and Sw are modified JONSWAP spectrums for 
% swell and wind peak, respectively.
% The energy is divided between the two peaks according
% to empirical parameters, which peak that is primary depends on parameters.   
% The empirical parameters are found for classes of Hm0 and Tp,
% originating from a dataset consisting of 20 000 spectra divided
% into 146 different classes of Hm0 and Tp. (Data measured at the
% Statfjord field in the North Sea in a period from 1980 to 1989.)
% The range of the measured  Hm0 and Tp for the dataset
% are from 0.5 to 11 meters and from 3.5 to 19 sec, respectively.
% See Torsethaugen (1996).
%
% Preliminary comparisons with spectra from other areas indicate that 
% some of the empirical parameters are dependent on geographical location.
% Thus the model must be used with care for other areas than the
% North Sea and sea states outside the area where measured data 
% are available.
%
% Example: [S,Sw,Ss]=torsethaugen([],[6 8],1)
%
% See also  jonswap, pmspec.

% References 
%  Torsethaugen, K. (1996)
%  Model for a doubly peaked wave spectrum 
%  Report No. STF22 A96204. SINTEF Civil and Environm. Engineering, Trondheim
%
%  Torsethaugen, K. (1994)
%  'Model for a doubly peaked spectrum. Lifetime and fatigue strength
%  estimation implications.' 
%  International Workshop on Floating Structures in Coastal zone,
%  Hiroshima, November 1994.
%
%  Torsethaugen, K. (1993)
%  'A two peak wave spectral model.'
%  In proceedings OMAE, Glasgow


% Tested on Matlab 5.3 
% History: 
% Revised pab april 2005
% revised pab 12.12.2000
% - fixed a bug: Hpw = 0 or Hps = 0 is now handled properly
% - sdata may now contain NaN's 
% revised pab 16.02.2000
%  -added ih, see at the end
% revised pab 20.12.1999
%  -added norm
%  -changed default value of Hm0 to 7 (equal to jonswap)
% revised pab 01.12.1999
%  - updated the reference
% revised  es 23.09.1999
%  - updated documentation
%  - added note
% by pab 23.06.1999

monitor=0;
% default values
Hm0   = 7;
Tp    = 11; 
Aj    = -1;
data2 = [Hm0, Tp, Aj];
nd2   = length(data2);
if nargin<3|isempty(plotflag),  plotflag = 0;       end
if (nargin>1) & ~isempty(sdata), 
  nd = length(sdata); 
  ind = find(~isnan(sdata(1:min(nd,nd2))));
  if any(ind), % replace default values with those from input data
    data2(ind) = sdata(ind);
  end
end
if (nd2>0) & (data2(1)>0), Hm0 = data2(1);end
if (nd2>1) & (data2(2)>0),  Tp = data2(2);end
if (nd2>2) & (data2(3)==0),  Aj = data2(3);end

w = [];
if nargin<1|isempty(w1), 
   wc = 33/Tp;
elseif length(w1)==1,
   wc = w1; 
else
   w = w1 ;
end
nw = 257;

if isempty(w),      w = linspace(0,wc,nw).'; end

w=w(:);

n      = length(w);
S      = createspec('freq');
S.S    = zeros(n,1);
S.w    = w;
S.norm = 0; % not normalized spectra
Sw = S;
Ss = S;

if Hm0>11| Hm0>max((Tp/3.6).^2, (Tp-2)*12/11)
  disp('Warning: Hm0 is outside the valid range')
  disp('The validity of the spectral density is questionable')
end
if Tp>20|Tp<3 
  disp('Warning: Tp is outside the valid range')
  disp('The validity of the spectral density is questionable')
end
wp = 2*pi/Tp;


%initializing

%sigma used in the peak enhancement factor
sa  = 0.07; % sigma_a for w/wp<1
sb  = 0.09; % sigma_b for w/wp>1

g = gravity; % m/s^2

% The parameter values below are found comparing the 
% model to average measured spectra for the Statfjord Field
% in the Northern North Sea.
Af  = 6.6;   %m^(-1/3)*sec
AL  = 2;     %sec/sqrt(m)
Au  = 25;    %sec
KG  = 35; 
KG0 = 3.5;
KG1 = 1;     % m
r   = 0.857; % 6/7
K0  = 0.5;   %1/sqrt(m)
K00 = 3.2;

M0  = 4;
B1  = 2;    %sec
B2  = 0.7;
B3  = 3.0;  %m
S0  = 0.08; %m^2*s
S1  = 3;    %m

% Preliminary comparisons with spectra from other areas indicate that 
% the parameters on the line below can be dependent on geographical location
A10 = 0.7; A1 = 0.5; A20 = 0.6; A2 = 0.3; A3 = 6;

Tf = Af*(Hm0)^(1/3);
Tl = AL*sqrt(Hm0);   % lower limit
Tu = Au;             % upper limit

%Non-dimensional scales
% New call pab April 2005
El = min(max((Tf-Tp)/(Tf-Tl),0),1); %wind sea
Eu = min(max((Tp-Tf)/(Tu-Tf),0),1); %Swell

%El = ((Tf-Tp)/(Tf-Tl)), %wind sea
%Eu = ((Tp-Tf)/(Tu-Tf)), %Swell


if Tp<Tf, % Wind dominated seas  
  % Primary peak (wind dominated)
  Nw  = K0*sqrt(Hm0)+K00;             % high frequency exponent
  Mw  = M0;                           % spectral width exponent
  Rpw = min((1-A10)*exp(-(El/A1)^2)+A10,1);
  Hpw = Rpw*Hm0;                      % significant waveheight wind
  Tpw = Tp;                           % primary peak period
  % peak enhancement factor
  gammaw = KG*(1+KG0*exp(-Hm0/KG1))*(2*pi/g*Rpw*Hm0/(Tp^2))^r;
  gammaw = max(gammaw,1);
  % Secondary peak (swell)
  Ns  = Nw;                % high frequency exponent
  Ms  = Mw;                % spectral width exponent
  Rps = sqrt(1-Rpw^2);
  Hps = Rps*Hm0;           % significant waveheight swell
  Tps = Tf+B1;
  gammas = 1;
  
 
  if Rps > 0.1
    disp('     Spectrum for Wind dominated sea')
  else 
    disp('     Spectrum for pure wind sea')
  end
else %swell dominated seas
   
  % Primary peak (swell)
  Ns  = K0*sqrt(Hm0)+K00;             %high frequency exponent
  Ms  = M0;                           %spectral width exponent
  Rps = min((1-A20)*exp(-(Eu/A2)^2)+A20,1);
  Hps = Rps*Hm0;                      % significant waveheight swell
  Tps = Tp;                           % primary peak period
  % peak enhancement factor
  gammas = KG*(1+KG0*exp(-Hm0/KG1))*(2*pi/g*Hm0/(Tf^2))^r*(1+A3*Eu);
  gammas = max(gammas,1); 
  
  % Secondary peak (wind)
  Nw   = Ns;                       % high frequency exponent
  Mw   = M0*(1-B2*exp(-Hm0/B3));   % spectral width exponent
  Rpw  = sqrt(1-Rps^2);
  Hpw  = Rpw*Hm0;                  % significant waveheight wind
  G0w  = Mw/((Nw/Mw)^(-(Nw-1)/Mw)*gamma((Nw-1)/Mw )); %normalizing factor
  if Hpw>0,
    Tpw  = (16*S0*(1-exp(-Hm0/S1))*(0.4)^Nw/( G0w*Hpw^2))^(-1/(Nw-1));
  else
    Tpw = inf;
  end
  %Tpw  = max(Tpw,2.5)
  gammaw = 1;
 
  if Rpw > 0.1
    disp('     Spectrum for swell dominated sea')
  else 
    disp('     Spectrum for pure swell sea')
  end    
end
if (3.6*sqrt(Hm0)<= Tp & Tp<=5*sqrt(Hm0))
   disp('     Jonswap range')
end
if monitor,
  disp(['Hm0 = ' num2str(Hm0)])
  disp(['Ns, Ms = ' num2str([Ns Ms]) '  Nw, Mw = ' num2str([Nw Mw])]) 
  disp(['gammas = ' num2str(gammas) ' gammaw = ' num2str(gammaw)])
  disp(['Rps = ' num2str(Rps) ' Rpw = ' num2str(Rpw)]) 
  disp(['Hps = ' num2str(Hps) ' Hpw = ' num2str(Hpw)])   
  disp(['Tps = ' num2str(Tps) ' Tpw = ' num2str(Tpw)]) 
end

%G0s=Ms/((Ns/Ms)^(-(Ns-1)/Ms)*gamma((Ns-1)/Ms )); %normalizing factor 

% Wind part

Sw = localJonswap(w,Hpw,Tpw,gammaw,Nw,Mw,Aj);
% Swell part
Ss = localJonswap(w,Hps,Tps,gammas,Ns,Ms,Aj);

 

S.S = Ss.S+Sw.S; % total spectrum
S.note=strcat('Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
if max(Ss.S)<max(Sw.S)
   Ss.note=strcat('Swell, secondary peak of Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
   Sw.note=strcat('Wind, primary peak of Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
else
   Ss.note=strcat('Swell, primary peak of Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
   Sw.note=strcat('Wind, secondary peak of Torsethaugen, Hm0=',num2str(Hm0),', Tp=',num2str(Tp));
end
if plotflag
  wspecplot(S,plotflag,'b')
  ih = ishold;
  hold on
  wspecplot(Ss,plotflag,'b--')
  wspecplot(Sw,plotflag,'b--')
  if ~ih, hold off, end
end

function  S = localJonswap(w,Hm0,Tp,gammai,N,M,A)
%local JONSWAP formulation

n1      = length(w);
S      = createspec('freq');
S.S    = zeros(n1,1);
S.w    = w;

if (Hm0<=0)
   return
end

%sigma used in the peak enhancement factor
sa  = 0.07; % sigma_a for wn<1
sb  = 0.09; % sigma_b for wn>1

wp  = 2*pi/Tp;
wn  = w/wp;

B  = (N-1)/M;
G0 = (N/M).^(B)*(M/gamma(B)); % Normalizing factor related to Pierson-Moskovitz form
G1 = (Hm0/4)^2/wp*G0;

% for w>wp
k = find(wn>1);
if any(k)
   S.S(k)=G1./(wn(k).^N).*(gammai.^(exp(-(wn(k)-1).^2 ...
                /(2*sb^2)))).*exp(-N/M*(wn(k)).^(-M));
end
% for 0<w<=wp
k = find(0<wn & wn<=1);
if any(k)
   S.S(k)=G1./(wn(k).^N).*(gammai.^(exp(-(wn(k)-1).^2 ...
      /(2*sa^2)))).*exp(-N/M*(wn(k)).^(-M));
end

if (gammai==1)
   A = 1;
end
if A<0, % normalizing by integration
   area = simpson(w,S.S);
   if area>0
      A = (Hm0/4)^2./area;% make sure m0=Hm0^2/16=int S(w)dw
   else
      A = 0;
   end
elseif A==0,% original normalization
   % NOTE: that  Hm0^2/16 generally is not equal to intS(w)dw
   %       with this definition of A if sa or sb are changed from the
   %       default values
   if 3<=N & N<=50 & 2<=M & M <=9.5 & 1<= gammai & gammai<=20
     f1NM = 4.1*(N-2*M^0.28+5.3)^(-1.45*M^0.1+0.96);
     f2NM = (2.2*M^(-3.3) + 0.57)*N^(-0.58*M^0.37+0.53)-1.04*M^(-1.9)+0.94;
     A = (1+ f1NM*log(gammai)^f2NM)/gammai;
   elseif N == 5 & M == 4,
     A = (1-0.287*log(gammai));
   else
      error('Not knowing the normalization because N, M or peakedness parameter was out of bounds!')
   end
end

S.S = S.S*A;


