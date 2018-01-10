function S = oscspec(sdata,z,w,s)
% OSCSPEC Spectral density for a harmonic oscillator
%         driven by white noise
%  
% CALL:  S = oscspec(sdata,z,w,s);
%
%        S     = the spectral density (structure array)
% 	 sdata = the data vector [wl wu n], where 
% 
%	    wl = lower truncation frequency  (default 4/257)
%	    wu = upper truncation frequency  (default 4)
%	    n  = number of evaluation points (default 257)
%
%       z,w,s = parameters in the equation (eq. 1) for the oscillator.
%               (default z=0.1, w=1, s=1)
%        
% Let W  be white noise, then the oscillator X is defined by
%  
%     X''(t) + 2wz X'(t) + w^2 X(t) = s W(t)          (eq. 1)
%   
% Example: 
%  data=[0.01 4 275];
%  S=oscspec(data,[],2.5);  % Peak frequency at w=2.5

% Tested on: Matlab 5.3
% History: 
% Correction by PJ 07-Jul-2005
%   Changed 'break' to 'return'
% revised es 25.05 00 small modifications of help text   
% Modified by jr 14.01.2000
% - updated check of nargins
% Modified by jr 12.01.2000
% - new names of variables and parameters
% - check of nargins introduced
% - updated documentation
% Modified by ir 11.01.2000
% - structure array introduced
% - parameters w,s allowed as input
% By Mats Frendahl 1993
  
if nargin<1 | isempty(sdata), sdata=[4/257 4 257];end
if nargin<2 | isempty(z), z=0.1;end
if nargin<3 | isempty(w), w=1;end
if nargin<4 | isempty(s), s=1;end

if (z<=0)|(z>=1)
  disp('   The parameter  z  must be in (0,1). Program will terminate.')
  return
end

wl=sdata(1); wu=sdata(2); n=sdata(3); 

wv=linspace(0,wu,n);

spv=s^2./((w^2-wv.^2).^2+(2*z*w)^2*wv.^2)/2/pi;

S=createspec;
S.S=spv; 
S.w=wv;
S.type='freq';
S.note='Spectrum, harmonic oscillator';
S.S(wv<wl)=0;
S.S(1)=0; % must be zero at zero freq since discrete spectrum
%S=floor(S*1e5+.5)/1e5;


