[filename, pathname] = uigetfile('*.dm3', 'Pick a DM3 file');
if ~(isequal(filename,0) || isequal(pathname,0))
	fn = fullfile(pathname,filename);
	[img sx sy]=ReadDM3_Matlab(fn);
	fprintf('read image from file %s\n',fn)
	fn2 = [fn(1:end-3),'img'];
	[filename2, pathname2] = uiputfile('*.img', 'Select a .img file',fn2);
    if ~(isequal(filename,0) || isequal(pathname,0))
		fn2 = fullfile(pathname2,filename2);
		binwrite2D(img.',fn2,10*sx,10*sy,0)
	end
end