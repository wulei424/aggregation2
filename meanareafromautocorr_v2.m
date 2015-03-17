function [A, radial_average_profile, radii_list, Cnorm, X, Y, dX, dY]=meanareafromautocorr_v2(I,maxr,X,Y,dX,dY)
% On first call 
% usage: [A, radial_average_profile, radii_list, Cnorm, X, Y, dX, dY]=meanareafromautocorr_v2(I,maxr);
% On any call afterwards
% usage: [A, radial_average_profile, radii_list, Cnorm]=meanareafromautocorr_v2(I,maxr,X, Y, dX, dY);
  
if nargin<2
maxr = 50; % round(0.4 * size(I,1));
end

% compute autocorrelation
I=double(I)-mean(double(I(:)));
R=xcorr2_fast(I);

dr=1;
center = (size(R)+3)/2;

if nargin ==2
    [out,X,Y,dX,dY]=radialaverageimage(R,center,dr,maxr);

else

[out,X,Y,dX,dY]=radialaverageimage(R,center,dr,maxr,X,Y,dX,dY);
end

radial_average_profile=out(:,1);
radii_list=out(:,2);

Cnorm = max(R(:));
r=find(radial_average_profile./Cnorm < 1/exp(1), 1, 'First');

if isempty(r)
    warning('correlation did not drop below threshold within sampling window. make maxlag larger or get larger images');
    A= pi* radii_list(end)^2;
else
    A= pi* radii_list(r)^2;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS


function [out,X,Y,dX2,dY2]=radialaverageimage(img,center,dr,maxr,X,Y,dX,dY,dostd, method)
    if nargin<9
        method=0;
    end
    
if nargin<3
    dr=1;
end

if nargin<4
    maxr=round(max(size(img)));
end

if nargin<8
    dostd=false;
end

xc=center(1);
yc=center(2);

if method==0
radii=dr:dr:maxr;

if nargin<7
    
    % Create the meshgrid to be used in resampling
[X,Y] = meshgrid(1:size(img,1),1:size(img,2));

for r=1:length(radii)
    radius=radii(r);
    % To avoid redundancy, sample at roughly 1 px distances
    num_pxls=2*pi*radius;
    theta=0:1/num_pxls:2*pi;
    dX{r}=radius*cos(theta);
    dY{r}=radius*sin(theta);
    %dX(:,r)=radius*cos(theta(:));
    %dY(:,r)=radius*sin(theta(:));
end
end


% radial slices via cellfun
%radial_slices=cellfun(@(x,y)interp2(X,Y,img,x+xc,y+yc),dX,dY,'UniformOutput',false);

% radial slices via arrayfun
if nargin<7
    dX2=cell2matNaN(dX);
    dY2=cell2matNaN(dY);
else
    dX2=dX;
    dY2=dY;
end
%radial_slices = arrayfun(@(r) interp2(X,Y,img,dX2(r,~isnan(dX2(r,:)))+xc, dY2(r,~isnan(dX2(r,:)))+yc), 1:length(radii),'UniformOutput',false);

% for cpu
% radial_average_profile=(cellfun(@nanmean, radial_slices));
% radial_sigma_profile=(cellfun(@nanstd, radial_slices));

% for gpuArray
% radial_average_profile=(cellfun(@mean, radial_slices,'UniformOutput',false));
% radial_sigma_profile=(cellfun(@std, radial_slices,'UniformOutput',false));
% radial_average_profile=gather([radial_average_profile{:}]);
% radial_sigma_profile=gather([radial_sigma_profile{:}]);

% radial average in one step
%dX2=cell2matNaN(dX);
%dY2=cell2matNaN(dY);

% interp2method
%interpmethod='*linear';
%interpfoo=@(r) nanmean(interp2q(X,Y,img,dX2(r,~isnan(dX2(r,:)))+xc, dY2(r,~isnan(dX2(r,:)))+yc, interpmethod));

%gridInterpolant method
interpmethod='linear';
F = griddedInterpolant(X',Y',img',interpmethod);
interpfoo=@(r) nanmean(F(dX2(r,~isnan(dX2(r,:)))+xc, dY2(r,~isnan(dX2(r,:)))+yc));

radial_slices = arrayfun(interpfoo, 1:length(radii),'UniformOutput',false);


radial_average_profile=gather([radial_slices{:}]);
radial_sigma_profile=[];

elseif method==1
%SLOW METHOD
% Create the meshgrid to be used in resampling
[X,Y] = meshgrid(1:size(img,1),1:size(img,2));
radii=dr:dr:maxr;
radial_sigma_profile=[];
for r = 1:length(radii)
% To avoid redundancy, sample at roughly 1 px distances
radius=radii(r);
num_pxls = 2*pi*radius;
theta = 0:1/num_pxls:2*pi;
x = xc + radius*cos(theta);
y = yc + radius*sin(theta);
sampled_radial_slice = interp2(X,Y,img,x,y);
radial_average_profile(r) = nanmean(sampled_radial_slice);

if dostd
    radial_sigma_profile(r) = nanstd(sampled_radial_slice);   
end

if nargin<7
    dX2=[];
    dY2=[];
else
    dX2=dX;
    dY2=dY;
end

end



end


out=[radial_average_profile(:), radial_sigma_profile(:), radii(:)];


function rmat=cell2matNaN(tcell)
% tcell = {[1,2,3], [1,2,3,4,5], [1,2,3,4,5,6], [1], []};  %# Sample array
% rmat=cell2matNaN(tcell);
% rmat =
% 
%      1     2     3   NaN   NaN   NaN
%      1     2     3     4     5   NaN
%      1     2     3     4     5     6
%      1   NaN   NaN   NaN   NaN   NaN
%    NaN   NaN   NaN   NaN   NaN   NaN
maxSize = max(cellfun(@numel,tcell));    %# Get the maximum vector size
fcn = @(x) [x nan(1,maxSize-numel(x))];  %# Create an anonymous function
rmat = cellfun(fcn,tcell,'UniformOutput',false);  %# Pad each cell with NaNs
rmat = vertcat(rmat{:});                  %# Vertically concatenate cells






function [xcorrM xlags ylags] = xcorr2_fast(IMa, IMb,biasflag,maxlags)
% xcorr2_fast(IMa,IMb)
% Calculates cross-correlation faster than Matlab's xcorr2
% C = xcorr2_fast (A, B, biasflag, maxlags)
%	Compute the 2D cross-correlation of matrices A and B.
% C = xcorr2_fast(A)
%      Compute two-dimensional autocorrelation of matrix A.
% C = xcorr2_fast(A,B)
%      Compute two-dimensional autocorrelation of matrix A, using no bias
%      correction, and no maxlag parameter

% C = xcorr2_fast(A,B, biasflag)
% where biasflag can be one of 'biased','unbiased','coeff',or 'none'
%      biased   - scales the raw cross-correlation by the maximum number
%                 of elements of A and B involved in the generation of 
%                 any element of C
%      unbiased - scales the raw correlation by dividing each element 
%                 in the cross-correlation matrix by the number of
%                 products A and B used to generate that element 
%      coeff    - normalizes the sequence so that the largest 
%                 cross-correlation element is identically 1.0.
%      none     - no scaling (this is the default).
%
% the default maxlags is none. When maxlag is passed, the correlation
% matrix is computed in full, but then truncated at the correct maxlag and
% so are the xlags and ylags truncated. This is provided as a convenience
% feature.
%
% 
% 8/15/2013
% updated biasflag options
% tests: figure;imagesc(xcorr2(I)-xcorr2_fast(I)) returns
% uncorrelated noise. You can also test by comparing against the 1D
% case by using a symmetric dataset.
% Carlos

% parse inputs
if (nargin < 1 || nargin > 4)
    disp('usage:"c = xcorr2(A [, B] [, scale])"');
end
if nargin == 1
    IMb = IMa;
    biasflag = 'none';
    maxlags=0;
elseif nargin == 2
    maxlags=0;
    if ischar(IMb)
        biasflag = IMb;
        IMb = IMa;
    else
        biasflag = 'none';
    end
elseif nargin==3
    maxlags=0;
end


%disp('using xcorr2_fast');
[Ha,Wa] = size(IMa);
[Hb,Wb] = size(IMb);


    Hc = Ha+Hb-1;
    Wc = Wa+Wb-1;
    

    Yc=(Hc+1)/2; %requires image to have even Hc and Wc
    Xc=(Wc+1)/2;
    

xcorrM = (real(ifft2(fft2(IMa,Hc,Wc).*fft2(rot90(IMb,2),Hc,Wc))));


if sum(maxlags)>0    
if length(maxlags)==2
    Hc = min(Hc-Yc,ensureodd(maxlags(1)));
    Wc = min(Wc-Xc,ensureodd(maxlags(2)));
elseif length(maxlags)==1
    Hc = min(Hc-Yc,ensureodd(maxlags));
    Wc = min(Wc-Xc,ensureodd(maxlags));
end
end

xcorrM=xcorrM((Yc-(Hc-1)/2):(Yc+(Hc-1)/2),(Xc-(Wc-1)/2):(Xc+(Wc-1)/2));
xlags=(-(Wc-1)/2):((Wc-1)/2);
ylags=(-(Hc-1)/2):((Hc-1)/2);

%xcorrM = fftshift(real(ifft2(fft2(IMa,Hc,Wc).*fft2(conj(IMb),Hc,Wc))));
% 
% FFT2(X) returns the two-dimensional Fourier transform of matrix X.
%     If X is a vector, the result will have the same orientation.
%  
%     FFT2(X,MROWS,NCOLS) pads matrix X with zeros to size MROWS-by-NCOLS
%     before transforming.


  %% bias routines by Dave Cogdell (cogdelld@asme.org)
  %% optimized by Paul Kienzle (pkienzle@kienzle.powernet.co.uk)
  %% optimized further by Carlos Ortiz (caortiz.phy@gmail.com)
  % see http://users.powernet.co.uk/kienzle/octave/matcompat/scripts/signal/xcorr2.m
  % and see http://www.mathworks.com/matlabcentral/fileexchange/78-xcorr2x-m
  if strcmpi(biasflag, 'biased'),
      xcorrM = xcorrM / ( min ([Ha, Hb]) * min ([Wa, Wb]) );
  elseif strcmpi(biasflag, 'unbiased'),
      
      lo = min ([Wa,Wb]); hi = max ([Wa, Wb]);
      row = [ 1:(lo-1), lo*ones(1,hi-lo+1), (lo-1):-1:1 ];
      lo = min ([Ha,Hb]); hi = max ([Ha, Hb]);
      col = [ 1:(lo-1), lo*ones(1,hi-lo+1), (lo-1):-1:1 ]';
      
      bias = col*row;
      mid=(size(bias)+1)/2;
      nlags=([Hc,Wc]-1)/2;
      xcorrM = xcorrM./bias((mid(1)-nlags(1)):(mid(1)+nlags(1)),(mid(2)-nlags(2)):(mid(2)+nlags(2)));
  elseif strcmpi(biasflag,'coeff'),
      xcorrM = xcorrM/max(xcorrM(:))';
  end
  
    function x=ensureodd(x)
        if mod(x,2)==0; 
             x=x+1;
        end
        
        
function Vq = interp2q(varargin)
%INTERP2 2-D interpolation (table lookup).
%
%   Some features of INTERP2 will be removed in a future release.
%   See the R2012a release notes for details.
%
%   Vq = INTERP2(X,Y,V,Xq,Yq) interpolates to find Vq, the values of the
%   underlying 2-D function V at the query points in matrices Xq and Yq.
%   Matrices X and Y specify the points at which the data V is given.
%
%   Xq can be a row vector, in which case it specifies a matrix with
%   constant columns. Similarly, Yq can be a column vector and it
%   specifies a matrix with constant rows.
%
%   Vq = INTERP2(V,Xq,Yq) assumes X=1:N and Y=1:M where [M,N]=SIZE(V).
%
%   Vq = INTERP2(V,K) returns the interpolated values on a refined grid 
%   formed by repeatedly dividing the intervals K times in each dimension.
%
%   Vq = INTERP2(V) is the same as INTERP2(V,1).
%
%   Vq = INTERP2(...,METHOD) specifies alternate methods.  The default
%   is linear interpolation.  Available methods are:
%
%     'nearest' - nearest neighbor interpolation
%     'linear'  - bilinear interpolation
%     'spline'  - spline interpolation
%     'cubic'   - bicubic interpolation as long as the data is
%                 uniformly spaced, otherwise the same as 'spline'
%
%   For faster interpolation when X and Y are equally spaced and monotonic,
%   use the syntax Vq = INTERP2(...,*METHOD).
%
%   Vq = INTERP2(...,METHOD,EXTRAPVAL) specificies a method and a scalar
%   value for Vq outside of the domain created by X and Y.  Thus, Vq will
%   equal EXTRAPVAL for any value of Yq or Xq which is not spanned by Y
%   or X respectively. A method must be specified for EXTRAPVAL to be used,
%   the default method is 'linear'.
%
%   All the interpolation methods require that X and Y be monotonic and
%   plaid (as if they were created using MESHGRID).  If you provide two
%   monotonic vectors, interp2 changes them to a plaid internally.
%   X and Y can be non-uniformly spaced.
%
%   For example, to generate a coarse approximation of PEAKS and
%   interpolate over a finer mesh:
%       [X,Y,V] = peaks(10); [Xq,Yq] = meshgrid(-3:.1:3,-3:.1:3);
%       Vq = interp2(X,Y,V,Xq,Yq); mesh(Xq,Yq,Vq)
%
%   Class support for inputs X, Y, V, Xq, Yq:
%      float: double, single
%
%   See also INTERP1, INTERP3, INTERPN, MESHGRID, scatteredInterpolant.

%   Copyright 1984-2013 The MathWorks, Inc.

narginchk(1,7); % allowing for an ExtrapVal

if (nargin == 7 || (nargin == 5 && ischar(varargin{end-1})) ) && ...
   (~isnumeric(varargin{end}) || ~isscalar(varargin{end}))
    error(message('MATLAB:interp2:InvalidExtrapval'))
end

% Parse the method and extrap val
[narg, method, ExtrapVal] = methodandextrapval(varargin{:});

% Construct the interpolant, narg should be 2,3 or 5 at this point.
if narg ~= 1 && narg ~= 2 && narg ~= 3 && narg ~= 5
    error(message('MATLAB:interp2:nargin'));
end
isspline = strcmpi(method,'spline');
if isspline 
    extrap = 'spline';
else
    extrap = 'none';
end
origvtype = 'double';
if narg == 1 || narg == 2
    %  interp2(V,NTIMES)
    V = varargin{1};
    [nrows,ncols] = size(V);
    if narg == 1
        ntimes = 1;
    else
        ntimes = floor(varargin{2}(1));
    end
    
    if ~isscalar(ntimes) || ~isreal(ntimes)
        error(message('MATLAB:interp2:NtimesInvalid'));
    end
    
    Xq = 1:1/(2^ntimes):ncols;
    Yq = (1:1/(2^ntimes):nrows)';
    [X, Y] = meshgridvectors(V);
    V = V.';
    if isspline
        [X, Y, V] = stripnanwrapper(X,Y,V); 
    end
    [V, origvtype] = convertv(V,method,Xq,Yq);
    F = makegriddedinterp({X,Y},V,method,extrap);
elseif narg == 3
    %  interp2(V,Xq,Yq)
    V = varargin{1};
    Xq = varargin{2};
    Yq = varargin{3};
    [X, Y] = meshgridvectors(V);
    V = V.';
    if isspline
        [X, Y, V] = stripnanwrapper(X,Y,V);
    end
    [V, origvtype] = convertv(V,method,Xq,Yq);
    F = makegriddedinterp({X,Y},V,method,extrap);
elseif narg == 5
    %  interp2(X,Y,V,Xq,Yq)
    X = varargin{1};
    Y = varargin{2};
    V = varargin{3};
    
    Xq = varargin{4};
    Yq = varargin{5};
    if isvector(X) && isvector(Y)
        V = V.';
        [X, Y, V] = checkmonotonic(X,Y,V);
        [V, origvtype] = convertv(V,method, Xq, Yq);
        if isspline
            [X, Y, V] = stripnanwrapper(X,Y,V);
        end
        F = makegriddedinterp({X, Y}, V, method,extrap);
    else
        V = V.';
        %[X, Y, V] = checkmonotonic(X',Y',V);
        X=X';
        Y=Y';
     %   V=V;
        [V, origvtype] = convertv(V,method, Xq, Yq);
        if isspline
            [X, Y, V] = stripnanwrapper(X,Y,V);          
        end
        F = makegriddedinterp(X, Y, V, method,extrap);
    end
end

if strcmpi(method,'cubic') && strcmpi(F.Method,'spline')
    % Uniformity condition not met
    gv = F.GridVectors;
    FV = F.Values;
    [X, Y, V] = stripnanwrapper(gv{:},FV);
    F = makegriddedinterp({X,Y},V,'spline');
end

% Now interpolate
iscompact = compactgridformat(Xq, Yq);
if iscompact || (isscalar(Xq) && isscalar(Yq) && strcmpi(F.Method,'spline'))
    Vq = F({Xq,Yq});
elseif isMeshGrid(Xq,Yq)
    iscompact = true;
    Xq = Xq(1,:);
    Yq = Yq(:,1);
    Vq = F({Xq,Yq});  
else
    Vq = F(Xq,Yq);
end

% impose the extrapolation value to the queries
% that lie outside the domain.
if ~isempty(ExtrapVal)
    Vq = imposeextrapval({Xq, Yq}, F.GridVectors, Vq, ExtrapVal, iscompact);
end

if  iscompact
    % Compact grid evaluation produces a NDGRID
    % Convert to MESHGRID
    Vq = Vq.';
end

if ~strcmp(origvtype,'double') && ~strcmp(origvtype,'single')
    Vq = cast(Vq,origvtype);
end
% end of interp2

% function stripnanwrapper
function [X, Y, V] = stripnanwrapper(X,Y,V)
[inan, jnan] = find(isnan(V));
ncolnan = length(unique(jnan));
nrownan = length(unique(inan));
if ncolnan == 0 && nrownan == 0
    return;
end
% Minimize loss of data. Strip rows instead of cols if there are less rows
if ncolnan > nrownan
    if isvector(X) && isvector(Y)
        % The X, Y and V are in compact NDGRID format
        [Y, X, V] = stripnansforspline(Y, X, V.'); % Swap on the way in & out
        V = V.';
    else
        [X, Y, V] = stripnansforspline(X',Y',V.');
        X = X.';
        Y = Y.';
        V = V.';
    end
else
    [X, Y, V] = stripnansforspline(X,Y,V);
end
warning(message('MATLAB:interp2:NaNstrip'));
if isempty(V) || isvector(V)
    error(message('MATLAB:interp2:NotEnoughPointsNanStrip'));
end

%function isMeshGrid
function isMG = isMeshGrid(X,Y)
if isvector(X) || ~ismatrix(X) || isempty(X) || ~isequal(size(X),size(Y))
    isMG = false;
elseif Y(1) ~= Y(1,end) || X(1) ~= X(end,1)  %quick check
    isMG = false;    
else
   isMG = norm(diff(X,[],1),1) == 0 && norm(diff(Y,[],2),1) == 0; 
end

%function convertv
function [V, origvtype] = convertv(V,method, Xq, Yq)
origvtype = class(V);
if ~isfloat(V) && (strcmpi(method,'nearest') || (isscalar(Xq) && isscalar(Yq)) )
    V = double(V);
end

%function makegriddedinterp
function F = makegriddedinterp(varargin)
try
    F = griddedInterpolant(varargin{:});
catch gime
    if iscell(varargin{1})
        method = varargin{3};
    else
        method = varargin{4};
    end
    if any(strcmp(gime.identifier,{'MATLAB:griddedInterpolant:NotNdGridErrId', ...
        'MATLAB:griddedInterpolant:NdgridNotMeshgrid2DErrId'}))
        error(message('MATLAB:interp2:InvalidMeshgrid'));
    elseif(strcmp(gime.identifier,'MATLAB:griddedInterpolant:DegenerateGridErrId') && strcmpi(method,'nearest'))
        error(message('MATLAB:interp2:DegenerateGrid'));
    else
        rethrow(gime);
    end
end

function [narg, method, ExtrapVal] = methodandextrapval(varargin)
%METHODANDEXTRAPVAL parses the method and ExtrapVal from the arguments
%   [NARG, METHOD, EXTRAPVAL] = METHODANDEXTRAPVAL(VARARGIN)
%   Parses VARARGIN to extract METHOD, and EXTRAPVAL and returns the number
%   of remaining arguments in NARG. This is a common utility function used
%   by INTERP2, INTERP3, INTERPN
%

%   Copyright 2012 The MathWorks, Inc.

narg = nargin;
ExtrapVal = [];
if isnumeric(varargin{end}) && isscalar(varargin{end}) && ischar(varargin{end-1})
    % User supplied an extrap val
    ExtrapVal = varargin{end};
    method_arg = varargin{end-1};
    narg = narg-2;
elseif ischar(varargin{end})
    method_arg = varargin{end};
    narg = narg-1;
else
    method_arg = 'linear';
end

% Parse the method_arg

if strncmp(method_arg,'*',1)
    method_arg = method_arg(2:end);
    % TODO issue a deprecation WARNING
end
if strncmpi(method_arg,'n',1)
    method = 'nearest';
elseif strncmpi(method_arg,'l',1) ||  (numel(method_arg) >=3 && strcmpi(method_arg(1:3),'bil'))
    method = 'linear';
elseif strncmpi(method_arg,'s',1)
    method = 'spline';
elseif strncmpi(method_arg,'c',1)  ||  (numel(method_arg) >=3 && strcmpi(method_arg(1:3),'bic') )
    method = 'cubic';
else
    method = method_arg;
end


function tf = compactgridformat(varargin)
% COMPACTGRIDFORMAT true if inputs are vectors of mixed orientation
%   TF = COMPACTGRIDFORMAT(X1, X2,X3,...Xn) returns true if X1, X2,X3,...Xn
%   are vectors and X1 and X2 have mixed orientations. This arrangement of 
%   vectors is used to implicitly define a grid in the INTERP2, INTERP3, and 
%   INTERPN functions.

%   Copyright 2012-2013 The MathWorks, Inc.

tf = all(cellfun(@isvector,varargin));
if tf && nargin > 1
    ns = cellfun(@(x)size(x)~=1, varargin, 'UniformOutput', false);
    tf = ~isequal(ns{:});
end


