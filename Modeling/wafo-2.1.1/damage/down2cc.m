function [cc,delta] = down2cc(mu,n)
%DOWN2CC Calculates the cycle count which has the highest damage 
%  given the downcrossing intensity. 
%
% CALL:  [cc,delta] = down2cc(cross,n);
%
%       cc    = the cycle count,
%       delta = the slice distance,
%       cross = a two column matrix with the levels  u  in the first
%               column and the number of downcrossings/downcrossing 
%               intensity in the second,
%       n     = the number of slice levels between 0 and maximum of
%               the number of downcrossings/downcrossing intensity.

%  Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

if exist('slice.dat'), delete slice.dat, end
if exist('crossint.dat'), delete crossint.dat, end

disp('   Writing data.')
fprintf('slice.dat','%4.0f\n',n)
fprintf('crossint.dat','%10.5f %10.5f\n',mu');

disp('   Starting Fortran executable.')
dos([ wafoexepath 'down2cc.exe']);

disp('   Loading data.')
load out.dat
cc=out(:,2:3);

delta=max(mu(:,2))/n;

delete slice.dat crossint.dat
