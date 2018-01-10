function options = rindoptset(varargin)
%RINDOPTSET Create or alter RIND OPTIONS structure.
%
%  CALL:  options = rindoptset(funcname,opts1,opts2,..,par1,val1,par2,val2,..);
%
%   options    = transformation options structure in which the named 
%                parameters have the specified values.  
%   funcname   = string giving the name of the function for which default
%                values for the options structure should be extracted.
%                Options are 'rind', 'spec2mmtpdf', 'spec2thpdf'.
%   opts1,
%   opts2..    = options structures
%   par1,par2..= strings identifying the parameter to alter
%   val1,val2..= corresponding values the parameters are altered to.
%   
%   RINDOPTSET combines the default options for a function given by FUNCNAME
%   with new options structures (OPTS1,OPTS2,...) and/or with the named
%   parameters (PAR1,PAR2,...) with the corresponding values (VAL1, VAL2,...).
%   The parameters are set in the same order as the input arguments.
%   Any parameters with non-empty values of the options struct overwrite
%   the corresponding old parameters. 
%   The input arguments can be given in any order with one exception:
%   PARx and VALx must be given in pairs in that order.
%   Any unspecified parameters for PARx are set to []. 
%   Parameters with value [] indicate to use the default value for that
%   parameter when OPTIONS is passed to the function. It is sufficient to
%   type only the 2 first characters to uniquely identify the parameter
%   or function name.  Upper case letters for parameter names and values
%   that are strings are ignored. If an invalid string is provided, the
%   default is used.
%   
%   RINDOPTSET with no input arguments and no output arguments displays all 
%   parameter names and their possible values.
%
%   RINDOPTSET with no input arguments creates an options structure
%   OPTIONS where all the fields are set to their default values.
%
%   
% RINDOPTSET PARAMETERS
%  METHOD  = INTEGER defining the integration method
%            0 Integrate by Gauss-Legendre quadrature  (Podgorski et al. 1999)
%            1 Integrate by SADAPT for Ndim<9 and by KRBVRC otherwise (default)
%            2 Integrate by SADAPT by Genz (1992) (Fast)
%            3 Integrate by KRBVRC by Genz (1993) (Fast)
%            4 Integrate by KROBOV by Genz (1992) (Fast)
%            5 Integrate by RCRUDE by Genz (1992)
%  XCSCALE = REAL to scale the conditinal probability density, i.e.,
%            f_{Xc} = exp(-0.5*Xc*inv(Sxc)*Xc + XcScale) (default XcScale =0)
%  ABSEPS  = REAL absolute error tolerance.       (default 0)
%  RELEPS  = REAL relative error tolerance.       (default 1e-3)
%  COVEPS  = REAL error tolerance in Cholesky factorization (default 1e-13)
%  MAXPTS  = INTEGER, maximum number of function values allowed. This 
%            parameter can be used to limit the time. A sensible 
%            strategy is to start with MAXPTS = 1000*N, and then
%            increase MAXPTS if ERROR is too large.    
%            (Only for METHOD~=0) (default 40000) 
%  MINPTS  = INTEGER, minimum number of function values allowed.
%            (Only for METHOD~=0) (default 0)
%  SEED    = INTEGER, seed to the random generator used in the integrations
%            (Only for METHOD~=0)(default floor(rand*1e9))
%  NIT     = INTEGER, maximum number of Xt variables to integrate
%            This parameter can be used to limit the time. 
%            If NIT is less than the rank of the covariance matrix,
%            the returned result is a upper bound for the true value
%            of the integral.  (default 1000)
%  XCUTOFF = REAL cut off value where the marginal normal
%            distribution is truncated. (Depends on requested
%            accuracy. A value between 4 and 5 is reasonable.)
%  XSPLIT  = parameters controlling performance of quadrature
%             integration:
%             if Hup>=xCutOff AND Hlo<-XSPLIT OR
%                Hup>=XSPLIT AND Hlo<=-xCutOff then
%             do a different integration to increase speed
%             in rind2 and rindnit. This give slightly different 
%            results
%            if XSPILT>=xCutOff => do the same integration allways
%            (Only for METHOD==0)(default XSPLIT = 1.5)   
%  QUADNO  = Quadrature formulae number used in integration of Xd
%            variables. This number implicitly determines number of nodes
%            used.  (Only for METHOD==0)
%  SPEED   = defines accuracy of calculations by choosing different 
%            parameters, possible values: 1,2...,9 (9 fastest,  default []).
%            If ~isempty(SPEED) the parameters, ABSEPS, RELEPS, COVEPS,
%            XCUTOFF, MAXPTS and QUADNO will be set according to INITOPTIONS.
%  
% Examples:
%  rindoptset('rind')
%
% See also  rind, initoptions

% History
% by pab 20.05.2003%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
% based on MATLAB's optimset


% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  help rindoptset
  return;
end

% Initialization
% Legal functions names
fnames = strvcat('rind','spec2mmtpdf','spec2thpdf','spec2tpdf2'); 
% Legal parameter names
names  = {'method','xcscale',...
	  'abseps','releps','coveps',...
	  'maxpts','minpts',...
	  'seed','nit','xcutoff',...
	  'xsplit','quadno', ...
	  'speed'};     
vals = {3,0,0,1e-3,1e-10,...
	40000,...
	0,...
	floor(rand*1e9),...
	1000,...
	[],...
	1.5,...
	[] ,...
	[]};

% Initialize options with default values
options = cell2struct(vals,names,2);
options = parseoptions(fnames,options,varargin{:});

if ~isempty(options.speed)
  options = initoptions(options.speed,options);
end
return
