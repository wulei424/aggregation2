function [A, radial_average_profile, radii_list, Cnorm, X, Y, dX, dY]=meanareafromautocorr_v2(I,maxr,X,Y,dX,dY)
% On first call 
% usage: [A, radial_average_profile, radii_list, Cnorm, X, Y, dX, dY]=meanareafromautocorr_v2(I,maxr);
% On any call afterwards
% usage: [A, radial_average_profile, radii_list, Cnorm]=meanareafromautocorr_v2(I,maxr,X, Y, dX, dY);
  
if nargin<2
maxr = 50; % round(0.4 * size(I,1));
end

% compute autocorrelation
R=xcorr2_fast(I);

dr=1;
center = (size(R)+3)/2;

[out,X,Y,dX,dY]=radialaverageimage(img,center,dr,X,Y,dX,dY);

radial_average_profile=out(:,1);
radii_list=out(:,2);

Cnorm = max(R(:));
r=find(radial_average./Cnorm < 1/exp(1), 1, 'First');

if isempty(r)
    warning('correlation did not drop below threshold within sampling window. make maxlag larger or get larger images');
    A= pi* radii_list(end)^2;
else
    A= pi* radii_list(r)^2;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS


function [out,X,Y,dX,dY]=radialaverageimage(img,center,dr,X,Y,dX,dY,dostd)
method=0;

if nargin<3
    dr=1;
end

if nagin<8
    dostd=false;
end

xc=center(1);
yc=center(2);

if method==0

if nargin<7
    
    % Create the meshgrid to be used in resampling
[X,Y] = meshgrid(1:size(img,1),1:size(img,2));
radii=dr:dr:round(max(size(img)));

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
radial_slices = arrayfun(@(r) nanmean(interp2(X,Y,img,dX2(r,~isnan(dX2(r,:)))+xc, dY2(r,~isnan(dX2(r,:)))+yc)), 1:length(radii),'UniformOutput',false);
radial_average_profile=gather([radial_slices{:}]);
radial_sigma_profile=[];

else
%SLOW METHOD
% Create the meshgrid to be used in resampling
[X,Y] = meshgrid(1:size(img,1),1:size(img,2));
radii=dr:dr:round(max(size(img)));
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