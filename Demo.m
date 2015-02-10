%% Source code for:
%% Article title: An automatic method for robust and fast cell detection in bright field images from high-throughput microscopy
%% MS ID        : 7277230959453875
%% Authors      : Felix Buggenthin, Carsten Marr, Michael Schwarzfischer, Philipp S Hoppe, Oliver Hilsenbeck, Timm Schroeder and Fabian J Theis
%% Journal      : BMC Bioinformatics, September 2013
%% When using this code in your publication, please cite accordingly
% 
%
% Copyright (C) 2013 Felix Buggenthin
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



%% define a parameter set
tiledim = 30;

lambda = 5; minSizeMSER = 30; maxSizeMSER = 4000; maxVariation = 1;

maxEcc = .7; minSizeSplit = 30; maxSizeSplit = 1000;


%load an image
I_org = imread('Demo1.png');

%run the code and visualize it
bw = segmentImage(I_org,'visualize',true);



