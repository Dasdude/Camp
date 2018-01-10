function y=binom(n,m)
% BINOM calculates the binomial coefficient n!/((n-m)!*m!)
%
% CALL:  y=binom(n,m);
%
% Example:%  
%    binom(5,2)    % Should be 10.

% tested on: matlab 5.x
% History:
% by pab 17.11.98
if 1,
  y=nchoosek(n,m);
else
  y=exp(gammaln(n+1)-gammaln(n-m+1)-gammaln(m+1));
end
