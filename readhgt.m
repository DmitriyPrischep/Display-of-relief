%READHGT import NASA SRTM data files (.HGT).
function varargout = readhgt(varargin)
sz3 = [1201,1201]; % SRTM3 tile size
novalue = intmin('int16'); % -32768
n = 1;

makeplot = 0;
decim = 0;
inter = 0;
nargs = 0;

if nargin == 0
	[filename,pathname] = uigetfile('*.hgt;*.hgt.zip','Select a HGT file');
	f = {[pathname,filename]};
	if filename == 0
		error('Please select a HGT file or use function arguments.');
	end
end

if nargin < (2 + nargs)
	lat = str2double(filename(2:3));
	if filename(1) == 'S'
		lat = -lat;
	end
	lon = str2double(filename(5:7));
	if filename(4) == 'W'
		lon = -lon;
	end

end

% pre-allocates X structure (for each file/tile)
X = repmat(struct('hgt',[],'lat',[],'lon',[]),[n,1]);

for n = 1:numel(f)
	% unzips HGT file if needed
	if ~isempty(strfind(f{n},'.zip'))
		X(n).hgt = char(unzip(f{n}));
		funzip = 1;
	else
		X(n).hgt = f{n};
		funzip = 0;
	end

	sz = sz3;
	if isempty(f{n})
		% offshore: makes a tile of zeros...
		X(n).z = zeros(sz);
	else
		% loads data from HGT file
		fid = fopen(X(n).hgt,'rb','ieee-be');
			X(n).z = fread(fid,'*int16');
		fclose(fid);
		switch numel(X(n).z)
		case prod(sz3)
			sz = sz3;
		otherwise
			error('"%s" seems not a regular SRTM data file or is corrupted.',X(n).hgt);
		end
		X(n).z = rot90(reshape(X(n).z,sz));

		% erases unzipped file if necessary
		if (funzip)
			delete(f{n});
		end
	end

	% builds latitude and longitude coordinates
	X(n).lon = linspace(lon(n),lon(n)+1,sz(2));
	X(n).lat = linspace(lat(n),lat(n)+1,sz(1))';
	
	% interpolates NaN (if not merged)
	if inter
		X(n).z = fillgap(X(n).lon,X(n).lat,X(n).z,novalue);
	end
end

if nargout == 0 || makeplot
    for n = 1:numel(X)
        fplot(X(n).lon,X(n).lat,X(n).z,decim,novalue)
    end
end

latitude = input('Input latitude: ');
longitude = input('Input longitude: ');

row = find(abs(X(n).lat - latitude) < 0.0001);
col = find(abs(X(n).lon - longitude) < 0.0001);

hgt = X(n).z(row, col);
disp(hgt)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fplot(x,y,z,decim,novalue)
%FPLOT plot the data using DEM function if exists, or IMAGESC

demoptions = {'latlon','legend','lake','nodecim'};

figure
if decim
	n = decim;
else
	n = ceil(sqrt(numel(z))/1201);
end
if n > 1
	x = x(1:n:end);
	y = y(1:n:end);
	z = z(1:n:end,1:n:end);
	fprintf('READHGT: In the figure data has been decimated by a factor of %d...\n',n);
end

if exist('dem','file')
	dem(x,y,z,demoptions{:})
else
	warning('For better results you might install the function dem.m')
	z(z==novalue) = 0;
	imagesc(x,y,z);
	if exist('landcolor','file')
		colormap(landcolor(256).^1.3)
	else
		colormap(jet)
	end
	% aspect ratio (lat/lon) is adjusted with mean latitude
	xyr = cos(mean(y)*pi/180);
	set(gca,'DataAspectRatio',[1,xyr,1])

	orient tall
	axis xy, axis tight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = fillgap(x,y,z,novalue)
%To optimize interpolation, we reduce the number of relevant data by building 
% a mask of surrounding pixels of novalue areas... playing with linear index!

sz = size(z);
k = find(z == novalue);
k(k == 1 | k == numel(z)) = []; % removes first and last index (if exist)
if ~isempty(k)
	[xx,yy] = meshgrid(x,y);
	mask = zeros(sz,'int8');
	k2 = ind90(sz,k); % k2 is linear index in the row order
	% sets to 1 every previous and next index, both in column and row order
	mask([k-1;k+1;ind90(fliplr(sz),[k2-1;k2+1])]) = 1; 
	mask(k) = 0; % removes the novalue index
	kb = find(mask); % keeps only border values
	z(k) = int16(griddata(xx(kb),yy(kb),double(z(kb)),xx(k),yy(k)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k2 = ind90(sz,k)

[i,j] = ind2sub(sz,k);
k2 = sub2ind(fliplr(sz),j,i); % switched i and j: k2 is linear index in row order
