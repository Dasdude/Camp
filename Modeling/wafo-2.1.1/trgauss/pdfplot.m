function H1=pdfplot(f,varargin)
%PDFPLOT Plot contents of pdf structures
%
% CALL: H = pdfplot(f,plotflag,x1,x2,x3,sym,method,shading)
%
%  H = handle to the created object
%
%  plot a PDF struct with the pollowing fields
%      
%      f.f      = pdf 
%      f.x      = cellarray of values in n dimensions (n=1:3) 
%  
% optional fields:   
%      f.labx   = cellarray of label strings          (n=1:3) 
%      f.title  = title string
%      f.cl     = contour levels for 2D PDF
%      f.pl     = Percent levels the given contour
%                 lines encloses.  
%     1D:
%      plotflag = 1 linear       (default)
%                 2 plot 1-f(x) on a semilog y-scale
%                 3 plot 1-f(x) on a log log scale
%                 11 plot F(x) = cumtrapz(x,f(x))
%                 12 plot 1-F(x) on a semilog y-scale
%                 13 plot 1-F(x) on a log log scale
%     2D:
%      plotflag = 1 contour plot (default)
%                 2 mesh
%                 3 surf
%                 4 waterfall
%                 5 pcolor
%                 6 contour3
%     3D: 
%      plotflag = 1 sliceomatic   (default)
%                 2 sliceomatic (with index x-,y-and z-labels)
%                 3 slice       
%                 4 contour f(X1,X2,X3(x1)), where x1 is an integer
%                 5 contour f(X1(x1),X2,X3), where x1 is an integer
%                 6 contour f(X1,X2(x1),X3), where x1 is an integer
%      x1,x2,x3 = are vectors defining where to slice for 3D data
%                 (default along the axis where f has its maximum)
%      sym      = plot symbol (default '-')  
%      method   = interpolation method for 3D slice 
%                 'linear' (default), 'cubic', or 'nearest'
%      shading  = controls the color shading of SURFACE and PATCH
%                 objects.  'faceted' (default), flat or 'interp'
%
%  Note: - sym,method and shading can be given anywhere after f and in
%          any order.
%
% See also   datastructures, qlevels, cltext

% Note: is only able to handle 1D,2D and 3D plot i.e. ndim=3

%Tested on: Matlab 5.3, 5.2
%History:
% revised pab March 2005
% -fixed some bugs
% revised pab 27.11.2002
% -added sliceomatic call -> reordered the plotflag order
% revised pab 07.01.2001
%  - fixed a bug for slice option 3D, reordered plotflag options for 1D.
% revised pab 23.11.2000
% - fixed a bug for calculation of contourlevels when f is not a pdf
% revised pab 08.11.2000
% - added the possibility that f is an array of structs
% - added plotflag 2,3 and 4 for 1D 
% revised pab 15.05.2000
%  - added shading, contour3
%  - now slicing along the axis where f has its maximum
% % revised pab 28.01.2000
%  - added point-and-click editing of all the text objects (title,
%    xlabel, ylabel) of the current figure
%  - improved the printing of contour level text and moved it into a
%    separate function cltext (this may be further improved)
%  - changed see also line
% revised es 24.01.2000 - dim=length(f.x) improving old, added som missing ;   
% revised pab 20.01.2000
%  - added pcolor
% revised pab 18.01.2000
%  - added return statement if plotflag==0
%  - added slice for 3D visualization
%  - added hold_state
%  - added pdfplotchk
% revised pab 5.11.1999
%  - changed PL to pl and CL to cl
% by pab 12.08.99

if ~isstruct(f) % secret option plot a matrix assuming the first column
  % is the independent variable
  [plotflag,sym] = pdfplotchk(varargin,1);
  switch plotflag
    case 0, return
    case 1,  H = plot(f(:,1),f(:,2:end),sym);
    case 2,  H = semilogy(f(:,1),1-f(:,2:end),sym);
    case 3,  H = loglog(f(:,1),-log(1-f(:,2:end)),sym);
    case 4,  H = semilogy(f(:,1),f(:,2:end),sym);
    case 5,  H = loglog(f(:,1),-log(f(:,2:end)),sym);
    case 6,  m = size(f,2)
             H = waterfall(f(:,1),1:m,f(:,2:end));
    case 11, H = plot(f(:,1),cumtrapz(f(:,1),f(:,2:end)),sym);  
    case 12, H = semilogy(f(:,1),1-cumtrapz(f(:,1),f(:,2:end)),sym);
    case 13, H = loglog(f(:,1),-log(1-cumtrapz(f(:,1),f(:,2:end))),sym); 
    otherwise, error('Unknown option for plotflag')
  end
  axis('square')
  set(gca,'FontSize',12)
  wafostamp;
   if (nargout>=1),     H1=H;   end
  return
end
hold_state = ishold; % remember old hold state
Nff=length(f);
if Nff>1
  cfig=gcf;
  for ix=1:Nff,
    if hold_state
      newplot
    else
      figure(cfig-1+ix)
    end
    pdfplot(f(ix),varargin{:})
  end
  return
end

cax  = newplot; % axes
cfig = get(cax,'Parent'); %cfig=gcf;

dim = length(f.x);
if dim>length(size(squeeze(f.f)))
  fsiz=size(f.f)
  dim=length(fsiz)-sum(fsiz==1); % Number of non-singleton dimensions
end
%dim,
[plotflag,sym,method,shad,x1,x2,x3] = pdfplotchk(varargin,dim,f);
if plotflag==0, return,end
switch dim			       
  case 1, %1D
    switch plotflag
      case 1,  H=plot(f.x{1},f.f,sym);
      case 2,  H=semilogy(f.x{1},1-f.f,sym);
      case 3,  H=loglog(f.x{1},-log(1-f.f),sym);
      case 11, H=plot(f.x{1},cumtrapz(f.x{1},f.f),sym);  
      case 12, H=semilogy(f.x{1},1-cumtrapz(f.x{1},f.f),sym);
      case 13, H=loglog(f.x{1},-log(1-cumtrapz(f.x{1},f.f)),sym);
    end
  case 2  %2D
    switch plotflag
      case {1,6,7,8,9},
	PL=0; 
	if isfield(f,'cl')&~isempty(f.cl) % check if contour levels is submitted
	  CL=f.cl;
	  if isfield(f,'pl'),PL=~isempty(f.pl);end % levels defines quantile levels? 0=no 1=yes
	else
	  CL=max(f.f(:))-range(f.f(:))*(1-[0.01 0.025 0.05 0.1 0.2 0.4 0.5 0.75]);
	  if 0 % automatic levels by using contours
	    c=contours(f.x{:},f.f.'); % calculate 8 levels
	    if isempty(c)
	      c=contours(f.x{:},f.f);%,7); % calculate levels
	    end
	    %CL = clevels(c);
	    limit = size(c,2);
	    ix = 1;
	    while(ix < limit)
	      CL(ix) = c(1,ix);
	      npoints = c(2,ix);
	      nexti = ix+npoints+1;
	      c(:,ix)=NaN;
	      ix = nexti;
	    end  
	    CL=unique(CL);
	  end
	end
	if PL,                            
	  clvec=sort(f.pl);                         
	else
	  clvec=sort(CL);
	end
	if any(plotflag==[1 8 9])
	  [cs hcs] = contour(f.x{:},f.f,CL,sym);
	else
	  [cs hcs] = contour3(f.x{:},f.f,CL,sym);
	end
	if any(plotflag==[1,6])
	  ncl=length(clvec);
	  if ncl>12, ncl=12; disp('   Only the first 12 levels will be listed in table.'),end
	  [hcl, axcl]=cltext(clvec(1:ncl),PL);  % print contour level text
	elseif any(plotflag==[7 9]) 
	  clabel(cs);
	else
	  clabel(cs,hcs);
	end
	
      case 2,	mesh(f.x{:},f.f); % meshz
      case 3,	surf(f.x{:},f.f);  %shading interp % flat, faceted       % surfc
      case 4,	waterfall(f.x{:},f.f);
      case 5,   pcolor(f.x{:},f.f); %shading interp % flat, faceted
      case 10,
	 [cs,hcs]=contourf(f.x{:},f.f); clabel(cs,hcs); fcolorbar(cs);
      otherwise, error('unknown option for plotflag')
    end
    if any(plotflag==[2:5])
       shading(shad);
    end
  case 3, %3D
    switch plotflag
     case 1,
      sliceomatic(f.x{:},f.f)
     case 2,
      sliceomatic(f.f)
     case 3,
      [X,Y,Z]=meshgrid(f.x{:});
      %method='linear';%, 'cubic','spline', or 'nearest'
      slice(X,Y,Z,f.f,x1,x2,x3,method);
      shading(shad);
     case 4,
      x1 = round(x1);
      contour(f.x{1:2},f.f(:,:,x1));
      if isempty(f.labx{3}),f.labx{3}='x3';end
      f.title=[f.title f.labx{3} ' = ' num2str(f.x{3}(x1))];
      f.labx{3}=[];
     case 5,
      x1 = round(x1);
      contour(f.x{[2 3]},squeeze(f.f(:,x1,:)).');
      if isempty(f.labx{1}),f.labx{1}='x1';end
      f.title=[f.title f.labx{1} ' = ' num2str(f.x{1}(x1))];
      f.labx{1}=[];
      f.labx=f.labx([2 3 1]);
     case 6,
      x1 = round(x1);
      contour(f.x{[1 3]},squeeze(f.f(x1,:,:)).');
      if isempty(f.labx{1}),f.labx{2}='x2';end
      f.title=[f.title f.labx{2} ' = ' num2str(f.x{2}(x1))];
      f.labx{2}=[];
      f.labx=f.labx([1 3 2]);
    end
end
if isfield(f,'labx') & max(size(f.labx))>=1
  Nf=max(size(f.labx));
  if Nf>=3, zlabel(f.labx{3}), end
  if Nf>=2, ylabel(f.labx{2}), end
  xlabel(f.labx{1})
end
if isfield(f,'title')
  title(f.title)
end

if exist('axcl','var'),% & ~isempty(axcl)
  set(cfig,'currentaxes',axcl(1))
  axis('square')
  set(cfig,'currentaxes',cax(1))
end
axis('square')
set(gca,'FontSize',12)
wafostamp;

%  The following two commands install point-and-click editing of
%   all the text objects (title, xlabel, ylabel) of the current figure:

set(findall(gcf,'type','text'),'buttondownfcn','edtext')
set(gcf,'windowbuttondownfcn','edtext(''hide'')')


if ~hold_state, 
   hold off, 
   %set(cfig,'NextPlot', 'replace'); 
end % reset to old hold state

if (nargout>=1)
  H1=H;
end
return



function [plotflag,sym,method,shad,x1,x2,x3] = pdfplotchk(P,dim,f)
%pdfplotCHK Helper function for pdfplot.
%
% CALL  [plotflag,sym,method,shad,x1,x2,x3]=pdfplotchk(P,dim) 
%
%   P = the cell array P of input arguments (between 0 and 6 elements)


% initialize output to default values
plotflag = 1;
sym='k-'; % Black dots is default
method='linear'; % linear is default
shad ='faceted';
if dim==3
  % Old call 
  %x1=mean(f.x{1});x2=mean(f.x{2});x3=mean(f.x{3});  
  % New call slicing where the maximum value is located
  [fmax, ind] = max(f.f(:));
  [I2,I1,I3] = ind2sub(size(f.f),ind);
  x1=f.x{1}(I1);x2=f.x{2}(I2);x3=f.x{3}(I3);
else
  x1=[];x2=[];x3=[];
end

Np=length(P);
try
  strix = cellfun('isclass',P,'char');
catch
  strix=zeros(1,Np);
  for ix=1:Np, % finding symbol strings 
    strix(ix)=ischar(P{ix});
  end
end
k=find(strix);
if any(k) % remove strings
  Nk=length(k);
  if Nk>3
    warning('More than 3 strings are not allowed in ')
  end
  for ix=1:length(k)
    switch lower(P{k(ix)})
    case {'flat','faceted','interp'}  , shad   = P{k(ix)};
    case {'linear','cubic', 'nearest'}, method = P{k(ix)};
    otherwise % plotsymbol is given
      sym = P{k(ix)};
    end
  end
  Np=Np-length(k);
  P={P{find(~strix)}}; % remove strings from input
end


if (Np>0) & ~isempty(P{1})
  plotflag=P{1};
end
if dim==3
  switch plotflag
    case {1,2,3}, % do nothing
    case 4, x1 = I3;
    case 5, x1 = I1;
    case 6, x1 = I2;
  end
end


if (Np>1) & ~isempty(P{2})
  x1=P{2};
end
if (Np>2) & ~isempty(P{3})
  x2=P{3};
end
 
if (Np>3) & ~isempty(P{4})
  x3=P{4};
end
return