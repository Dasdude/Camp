function F = wfrechcdf(x,a,c)
%WFRECHCDF Frechet cumulative distribution function
%
% CALL:  F = wfrechcdf(x,a,c);
%
%        F = density function evaluated at x
%     a, c = parameters
%
% The Frechet distribution is defined by its cdf
%
%  F(x;a,c) = exp(-(x/a)^(-c)), x>=0, a,c>0
%
% Example: 
%   x = linspace(0,6,200);
%   F1 = wfrechcdf(x,1,1); F2 = wfrechcdf(x,2,2); F3 = wfrechcdf(x,2,5);
%   plot(x,p1,x,F2,x,F3)

% Reference: 

% Tested on; Matlab 5.3
% History: 
% Added PJ 10-May-2001


error(nargchk(3,3,nargin))

[errorcode, x, a, c] = comnsize (x,a, c);
if (errorcode > 0)
  error ('x, a and c must be of common size or scalar');
end

F=zeros(size(x));

ok = ((c > 0)  & (a > 0));

k = find (x>=0&ok);
if any (k)  
  F(k)=exp(-(x(k)./a(k)).^(-c(k)));
end

k1 = find (~ok);
if any (k1)
  tmp=NaN;
  F(k1) = tmp(ones(size(k1)));
end
