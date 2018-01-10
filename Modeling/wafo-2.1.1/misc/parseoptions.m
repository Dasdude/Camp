function options = parseoptions(varargin)
%PARSEOPTIONS Create or alter a OPTIONS structure.
%
%  CALL:  options = parseoptions(legalFnames,options,funcname,opts1,opts2,...,
%                                par1,val1,par2,val2,...);
%         options = parseoptions(options,opts1,opts2,...,
%                                par1,val1,par2,val2,...);  
%
%   options    = OUT: options structure in which the named 
%                     parameters have the specified values.  
%                IN:  options structure giving the legal names and
%                     default options  
%   legalFnames = character array giving the names of functions for
%                 which default values for the options structure could be
%                 extracted.
%   funcname   = string giving the name of the function for which default
%                values for the options structure should be extracted.
%                (Must be equal to one of the names in LEGALFNAMES.)  
%   opts1,
%   opts2..    = options structures
%   par1,par2..= strings identifying the parameter to alter
%   val1,val2..= corresponding values the parameters are altered to.
%   
%   PARSEOPTIONS combines the default options for a function given by
%   FUNCNAME with new options structures (OPTS1,OPTS2,...) and/or with
%   the named parameters (PAR1,PAR2,...) with the corresponding values
%   (VAL1, VAL2,...).
%   The parameters are set in the same order as the input arguments.
%   Any parameters with non-empty values of the options struct overwrite
%   the corresponding old parameters. 
%   The input arguments can be given in any order with one exception:
%   PARx and VALx must be given in pairs in that order.
%   Any unspecified parameters for PARx are set to []. 
%   Parameters with value [] indicate to use the default value for that
%   parameter when OPTIONS is passed to the function. It is sufficient to
%   type only the first characters that uniquely identify the parameter
%   or function name.  Upper case letters for parameter names and values
%   that are strings are ignored. If an invalid string is provided, the
%   default is used.
%   
%   PARSEOPTIONS with no input and no output arguments displays this help 
%
% Examples:
%  defaultoptions = struct('test',[],'integrate',[] )
%  parseoptions(defaultoptions,'int','yes')
%  opt = defaultoptions;
%  opt.test = 'yes';  
%  parseoptions(defaultoptions,'int','yes',opt)
  
% History
% by pab 14.05.2003
% based on MATLAB's optimset

if (nargin == 0) 
  % Display help
  help parseoptions
  return;
end

% Initialization
fnames = '';

% Legal functions names
ix = 1;
if ((ix<=nargin) & (ischar(varargin{ix}))),
  fnames = varargin{ix};
  ix = 2;
end
if ((ix<=nargin) & (isstruct(varargin{ix}))),
  options = varargin{ix};
  ix = ix + 1;
else
  error('Struct expected.')
end

% Legal parameter names
namesc = fieldnames(options);
names  = lower(strvcat(namesc{:}));

expectval = 0;         % start expecting a name or stucture, not a value
while ix <= nargin
  arg = varargin{ix};
  if expectval
    options   = setfield(options,namesc{iy},arg);
    expectval = 0;
  else
    switch lower(class(arg))
      case 'char',
	lowArg = strtok(lower(arg),'.');      % Downcase and strip .m extensions
	iy = findlname(lowArg,fnames,0);      % Find legal function name
	if length(iy)==1,
	  opt1    = defaultopt(fnames(iy,:)); % Get default options for function
	  options = combineopt(options,opt1); % Combine with old options
	else,                                
	  iy = findlname(lowArg,names,1);     % Find legal parameter name
	  if isempty(iy),                     % Illegal parameter name
	    ix = ix+1;                        % Skip next input  
	  else
	    expectval = 1;                    % expect a value next
	  end
	end  
      case 'struct'
	options = combineopt(options,arg);
      otherwise
	error(sprintf('Expected argument %d to be a string or a struct.', ix));
    end
  end
  ix=ix+1;
end

if expectval
   error(sprintf('Expected value for parameter ''%s''.', arg));
end
return

function iy = findlname(arg,names,chkempty)
% FINDLNAME find index to legal parameter name  

iy = strmatch(arg,names);
%iy = find(strncmp(lowArg,fnames,2));
Niy = length(iy);
if Niy<1 & chkempty                      % No matches
  msg = ['Unrecognized parameter name ' arg '.'];
  warning(msg);
elseif Niy > 1      % More than one match
  % Find exact matches in case any names are subsets of others
  k = strmatch(arg,names,'exact');
  if length(k) == 1
    iy = k;
  else
    msg = ['Ambiguous parameter name ' arg ' ('];
    for k = iy'
      msg = [msg ', ' deblank(names(k,:))];
    end
    msg = [ msg ').']
    error(msg);
  end
end
return

function options = combineopt(options,newopts)
%COMBINEOPT combines an existing options structure
%   OPTIONS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%   with non-empty values overwrite the corresponding old parameters in
%   OPTIONS. 

if ~isempty(newopts)                      % [] is a valid options argument
  newNames    = fieldnames(newopts);
  legalNames  = fieldnames(options);
  
  [commonNames, ia, ib] = intersect(lower(newNames),lower(legalNames));
  if isempty(ia),
    return;
  end
  cval = struct2cell(newopts);
  ind  = findNonEmptyCells(cval(ia));
  if any(ind)
    val  = struct2cell(options);
    val(ib(ind)) = cval(ia(ind));
    options = cell2struct(val,legalNames,1);
  end
end
return


function options = defaultopt(fname)
% DEFAULTOPT Get default OPTIONS structure from the function FNAME

fname = deblank(fname);

if ~exist(fname)
  msg = ['No default options available: the function ' fname ' does not exist on the path.'];
  error(msg)
end
try 
  options = feval(fname,'defaults');
catch
  msg = ['No default options available for the function ' fname '.'];
  error(msg)
end
return

function ind = findNonEmptyCells(carray)
  try, % matlab 5.3 or higher
    ind = find(~cellfun('isempty',carray)).';
  catch
    % Slow 
    n = length(carray);
    ind1 = zeros(1,n);
    for ix = 1:n
      ind1(ix) = isempty(carray{ix});
    end
    ind = find(~ind1);
  end
  return