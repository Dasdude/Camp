function ind=findrfc(tp,h)
%FINDRFC Finds indices to rainflow cycles of a sequence of TP.
%
% CALL:  RFC_ind = findrfc(TP,h);
%
%        TP  = vector of turningpoints (NB! Only values, not sampled times)
%
%    RFC_ind = indices to the rainflow cycles of the 
%		original sequence TP.
%
%         h  = a threshold, must be larger than zero. 
%              if h>0, then all rainflow cycles with height 
%              smaller than  h  are removed.
%
%  This function is not implemented as a matlab function; instead, a 
%  mex file (originally written in C) is utilized.
% 
%  Example:
%   x = load('sea.dat'); 
%   tp = dat2tp(x); 
%   ind = findrfc(tp(:,2),0.3); 
%   waveplot(x,tp(ind,:),1,1)  
% 
%  See also  tp2rfc, dat2tp. 
  
% This is a modified version of rfcfilt (found in WAT), which is about 20 
% to 30 times faster than rfcfilt (on a PentiumII 233 MHz  
% with 32 MB ram and Matlab 5.0 under Linux). The reason is
% that this version does not save TP to disk. Instead it passes 
% the arguments directly to the executeable file. 
% However, this solution requires different input and returns
% indices to the rfc instead of the rfc itself.
% It also ignores the first turningpoint if that is a maximum and 
% starts on the first minimum when finding the sequence of rfc. 

% Tested on Matlab 6.1,6.0, 5.2
% History:
% revised pab 10.08.2003
% - fixed a bug in the example  
% Revised by pab 24.07.1999


ind=[];
disp('FINDRFC is not implemented as a m-function')
disp('                   compile the mexfile findrfc.c before you try again.')
error('findrfc error')
return