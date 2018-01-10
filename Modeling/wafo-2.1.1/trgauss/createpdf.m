function f=createpdf(ndim)
% CREATEPDF PDF class constructor
%
% CALL:  f=createpdf(ndim)
%
%  creates a struct with the following fields
%      
%      f.f      = pdf (ndim-dimensional matrix)
%      f.x      = cellarray of vectors of lags in ndim dimensions (ndim cells)
%      f.labx   = cellarray of strings of label strings (ndim cells)           
%      f.title  = title string
%      f.note   = note string
%      f.date   = creation date and time
%
%      ndim     = # dimensions (default 1)
%
%  Examples:   f = createpdf(2) gives the structure
%        f: []
%        x: {2x1 cell}
%     labx: {2x1 cell}
%    title: []
%     note: []
%     date: '16-Oct-1999 16:56:54'
%
% See also  datastructures, pdfplot, createcov

%Tested on: Matlab 5.3
%History:
% revised es 25.10.1999 help and cosmetics 
% revised pab 16.10.1999
% changed .xi and .labxi to .x and .labx cellarrays for 
%      faster computation and easier access
% by pab 16.09.99


if nargin<1|isempty(ndim)
 ndim=1;
end

f=struct('f',[]);
f.x=cell(ndim,1);
f.labx=cell(ndim,1);
f.title=[];
f.note=[];
f.date=datestr(now);

