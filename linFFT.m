function [frames,frameFFT,nframes,sectorsum,peaks,peak_ft,peak_factor_ft,...
    backtransf_max,backtransf_ft,wavelength_max,wavelength_ft,map_max,map_ft] ... 
    = linFFT(image,framesize,threshold,dispfig,step,tmin,tmax,resolution)
% linFFT is a function for detecting and analysing linear features in an 
% image using the Fast Fourier Transform. 
%   
%   image       image to perform this function on. can be unedited.
%
%   framesize   size of sub-images (in px, should be either 51 or 101)
%               default = 101 px
%
%   resolution  spatial resolution of image, in m/px.
%               default = 1.06 m/px (resolution of Hathor image)
%
%   step        angular step of the radial sector, in degrees. 
%               default = 5 degrees
%               larger step = less precision in the iFFT images
%               smaller step = lower peaks, less certainty in peak detection
%
%   threshold   minimum height of a maximum in the intensity spectrum to 
%               count as a peak, relative to the mean intensity of all sectors
%               (e.g. "1.5" means 150% of the mean)
%               default = 1.5
%
%   tmin        min. theta angle of features that are considered
%               default = 141.5 degrees (works best on Hathor features)
%
%   tmax        max. theta angle of features that are considered
%               default = 166.5 degrees (works best on Hathor features)
%
%   dispfig     toggles whether function outputs are shown in a figure
%               = -1    displays two maps made from backtransformations of
%                       the sub-images (default)
%               =  0    no figure, no text. results are passed to output
%                       variables
%               =  0.1  no figure, results are shown as text and passed to
%                       output variables
%               >= 1    index of sub-image for which a multiplot is shown to display 
%                       the sub-image's results

%% set up image and parameters

% check and adjust image's dimensions to end in '01':
verticalDimMod = mod(size(image,1),100); 
horizontalDimMod = mod(size(image,2),100);
switch verticalDimMod
    case 0
        image(end+1,:,:) = 0;    
    case num2cell(2:79) % up to '79', crop image
        verticalDimNew = (size(image,1))-(verticalDimMod-1);
        image = image(1:verticalDimNew,:,:);
        fprintf('Note: image was cropped at the bottom by %i pixels.\n',...
            verticalDimMod-1);
    case num2cell(80:99) % > '80' add black padding
        verticalDimNew = (size(image,1))+(101-verticalDimMod);
        image(end:verticalDimNew,:,:) = 0;
end
switch horizontalDimMod
    case 0
        image(:,end+1,:) = 0;    
    case num2cell(2:79)
        horizontalDimNew = (size(image,2))-(horizontalDimMod-1);
        image = image(:,1:horizontalDimNew,:);
        fprintf('Note: image was cropped on the right side by %i pixels.\n',...
            horizontalDimMod-1);
    case num2cell(80:99)
        horizontalDimNew = (size(image,2))+(101-horizontalDimMod);
        image(:,end:horizontalDimNew,:) = 0;
end
originalimage = image; % saves a copy of the original image (in color) for display
image = im2double(image);
if size(image,3) > 1 % if image is NOT greyscale          
    image = rgb2gray(image);
end

if ~exist('resolution','var') 
    resolution = 1.06; % spatial resolution of image (in m/px)
end
if ~exist('step','var')
    step = 5; % angular step of the radial sector (in degrees)
end
if ~exist('framesize','var')
    framesize = 101; % size of analysed sub-images (in px)
end
    fz = framesize;
if ~exist('threshold','var')
    threshold = 1.5; % min. intensity for peak-detection of a sector,
                     % relative to the mean intensity of all sectors
end
if ~exist('tmin','var')
    tmin = 141.5; % min. angle (theta) of target features
end
if ~exist('tmax','var')
    tmax = 166.5; % max. angle (theta) of target features
end
if ~exist('dispfig','var')
    dispfig = -1; % by default, two maps are generated
end

rmin = 6;                               % minimum radius
rmax = round(0.4 * fz);                 % maximal radius
nsec = 360/step;                        % number of sectors
sectors = [0:(step/2):360 (step/2)]';   % vector of overlapping sector boundaries, 
                                        % ends on 'step/2'
c = round(fz/2); % centre pixel of frame

window = mat2gray(fspecial('Gaussian',fz,0.4*fz)); % mask to reduce border effects

% for each pixel of a (fz x fz) frame, find theta-value (t) and radius (r):
t = zeros(fz); r = zeros(fz);
for i = 1:fz
    for j = 1:fz
        [th,ra] = cart2pol(-(i-c),j-c);
        if rad2deg(th) >= 0
            t(i,j) = rad2deg(th);
        else
            t(i,j) = rad2deg(th) + 360;
        end
        r(i,j) = ra;
    end
end
clear('ra','th','i','j')

% identify which pixels belong to each sector (pre-define this for speed):
sectormasks = cell(2*nsec,1); [sectorpixelsum,sps] = deal(zeros(2*nsec,1));
for k = 1:(2*nsec)-1
    m = zeros(fz);
    for j = 1:fz
        for i = 1:fz
            if (t(i,j) >= sectors(k)) && (t(i,j) < sectors(k+2)) ...
                    && (r(i,j) > rmin) && (r(i,j) <= rmax)
               m(i,j) = 1;
            end
        end
    end
    sectormasks{k} = m;
    sectorpixelsum(k) = sum(sum(sectormasks{k}));
    sps(k) = sectorpixelsum(k)/max(sectorpixelsum); % scaling factor to avoid 
                                                    % uneven #of px per sector
end
for k = 2*nsec      % special case: crossing zero angle (north)
    m = zeros(fz);
    for j = 1:fz
        for i = 1:fz
            if (t(i,j) >= sectors(k)) && (t(i,j) <= 360) ...
                    && (r(i,j) > rmin) && (r(i,j) <= rmax)
               m(i,j) = 1;
            end
            if (t(i,j) >= 0) && (t(i,j) < sectors(k+2)) ...
                    && (r(i,j) > rmin) && (r(i,j) <= rmax)
               m(i,j) = 1;
            end
        end
    end
    sectormasks{k} = m;
    sectorpixelsum(k) = sum(sum(sectormasks{k}));
    sps(k) = sectorpixelsum(k)/max(sectorpixelsum); 
end 
% note: sectormask{k} contains pixels between sectors(k) and sectors(k+2)
clear('m','sectorpixelsum','i','j')

%% define sub-images ('frames')

nframes_h = round((size(image,2))/fz) ;     % # of horizontal frames, non-overlapping
nframes_v = round((size(image,1))/fz) ;     % # of vertical frames, non-overlapping
nframes = (2*nframes_h-1) * (2*nframes_v-1);  % total number of frames, overlapping
frames = cell(nframes,fz);

f = 1;
for v = 1:0.5:nframes_v
    for h = 1:0.5:nframes_h
        updown = ((v-1)*(fz-1))+1:(v*(fz-1))+1;
        leftright = ((h-1)*(fz-1))+1:(h*(fz-1))+1;
        frames{f} = image(updown,leftright); 
        %(if code fails in line 158, the image has wrong dimensions)
        f = f+1;
    end
end
clear('f','updown','leftright','v','h')

%% compute FFT of all sub-images ('frames')

[frameFFT,frameFFT2] = deal(cell(nframes,1));                         
             
for i = 1:nframes
    frameFFT{i} = fftshift(fft2(frames{i}.*window));
    frameFFT2{i} = fftshift(fft2(frames{i}.*window)).^2;  % squared to enhance peaks
    frameFFT2{i}(c-2:c+2,c-2:c+2) = 0;  % blots out center pixels to enhance display
                                        % (does not affect calculations because of rmin)
end

%% find direction of FFT-intensity-peaks

% find cumulative FFT-intensities per radial sector:
sectorsum = zeros(nframes,2*nsec);  % (each row one frame, each column one sector)
[intensity,peaks,peakFactors,peakSectors] = deal(zeros(nframes,4));
[peak_ft_index,peak_ft,peak_factor_ft,peakSector_ft] = deal(cell(nframes,1));
 
for i = 1:nframes
    for k = 1:(2*nsec)
        sectorsum(i,k) = sum(sum(sectormasks{k}.*abs(frameFFT2{i})));
            % sum of FFT intensities for each sector
            % be careful with indices: sectorsum(k) contains sum of intensities for modes 
            % whose theta value lies between sectors(k) and sectors(k+2)
            % and whose radius lies between rmin and rmax
        sectorsum(i,k) = sectorsum(i,k) / sps(k); % normalize for #of pixels in sector
    end
 
    [ints,pks] = findpeaks(sectorsum(i,1:nsec+1),sectors(2:nsec+2),...
        'SortStr','descend',...
        'MinPeakDistance',15,'MinPeakWidth',5,...       
        'MinPeakHeight',threshold*mean(sectorsum(i,:))); 

    intensity(i,1:length(ints)) = ints; 
    peaks(i,1:length(pks)) = pks';
    
    for s = 1:length(ints) % determine which sectors the peaks belong to
        peakSectors(i,s) = find(sectors == peaks(i,s))-1;
    end
    clear('s')
    
    peak_ft_index{i} = find(peaks(i,:)>=tmin & peaks(i,:)<=tmax,1,'first'); 
        % these ft-related variables must be handled as cells to accommodate
        % empty entries if no peak is found!
    peak_ft{i} = peaks(i,peak_ft_index{i});
        % peaks are sorted by descending intensity (see above), so if more than
        % one peak_ft were to be found within [tmin,tmax], the first value is
        % more relevant (as it has higher intensity)
        
    peakSector_ft{i} = find(sectors == peak_ft{i})-1;
    
    clear('s')
    
    peakFactors(i,1) = intensity(i,1)/mean(sectorsum(i,:)); 
    peakFactors(i,2) = intensity(i,2)/mean(sectorsum(i,:)); 
    peakFactors(i,3) = intensity(i,3)/mean(sectorsum(i,:)); 
    peakFactors(i,4) = intensity(i,4)/mean(sectorsum(i,:)); 
    peak_factor_ft{i} = intensity(peak_ft_index{i})/mean(sectorsum(i,:)); 

end
clear('k','sps','ints','pks')

%% compute backtransformations

placeholder = zeros(fz); placeholder(2:end-1,2:end-1) = 1;
[fft_reduced,fft_reduced_ft,backtransf_max,backtransf_ft] = deal(cell(nframes,1));
keep = zeros(nframes,2);

% max peak:

for i = 1:nframes
    fft_reduced{i} = frameFFT{i};
    if peaks(i,1) ~= 0
        keep(i,:) = [peaks(i,1)-(step/2) peaks(i,1)+(step/2)];            
        for jj = 1:fz
            for ii = 1:fz
                if ((t(ii,jj) >= keep(i,1)) && (t(ii,jj) <= keep(i,2))) || ...
                        ((t(ii,jj) >= keep(i,1)+180) && (t(ii,jj) <= keep(i,2)+180))
                    % for all px inside of 'keep', do nothing
                else
                    fft_reduced{i}(ii,jj) = 0; % set all px outside of 'keep' to 0
                end
            end
        end
        backtransf_max{i} = ifft2(ifftshift(fft_reduced{i}));
        
        % alternative backtransformation based on sectormasks (excludes
        % Fourier modes outside of [rmin rmax]). This backtransformation
        % visualises which wavelengths are most relevant for finding \theta.
        fft_reduced_alt{i} = frameFFT{i} .* sectormasks{peakSectors(i,1)}; %#ok<AGROW>
        backtransf_max_alt{i} = real(ifft2(ifftshift(fft_reduced_alt{i}))); %#ok<AGROW>
        
    else
        backtransf_max{i} = placeholder;
    end
end

% feature peak: 
for i = 1:nframes
    fft_reduced_ft{i} = frameFFT{i};
    if ~isempty(peak_ft{i})
        keep(i,:) = [peak_ft{i}-(step/2) peak_ft{i}+(step/2)];            
        for jj = 1:fz
            for ii = 1:fz
                if ((t(ii,jj) >= keep(i,1)) && (t(ii,jj) <= keep(i,2))) || ...
                        ((t(ii,jj) >= keep(i,1)+180) && (t(ii,jj) <= keep(i,2)+180))
                    % for all px inside of 'keep', do nothing
                else
                    fft_reduced_ft{i}(ii,jj) = 0; % set all px outside of 'keep' to 0
                end
            end
        end
        backtransf_ft{i} = ifft2(ifftshift(fft_reduced_ft{i}));
    else
        backtransf_ft{i} = placeholder;
    end
end

clear('placeholder')

%% FFT power spectra in direction of peaks

[power_peak1,power_peak2,power_peak3,power_peak4,power_peak_ft] ...
    = deal(cell(nframes,1));
wavelength_max = zeros(nframes,4); wavelength_ft = zeros(nframes,1);
keep = zeros(nframes,2);

for i = 1:nframes
    if peaks(i,1) ~= 0
        keep(i,:) = [peaks(i,1)-(step/2) peaks(i,1)+(step/2)];  
        count = 1; % counts through pixels (needed as index)
        for jj = 1:fz         % horizontal
            for ii = 1:fz     % vertical
                if (t(ii,jj) >= keep(i,1) && (t(ii,jj) <= keep(i,2)) ...
                            && (r(ii,jj) >= rmin) && (r(ii,jj) <= rmax))
                % for each pixel within 'keep', determine intensity   
                    power_peak1{i}(count,1) = r(ii,jj);                    % radius of pixels
                    power_peak1{i}(count,2) = fz/power_peak1{i}(count,1);  % wavelength of Fourier mode of pixels
                    power_peak1{i}(count,3) = abs(frameFFT2{i}(ii,jj));    % power of Fourier mode
                    power_peak1{i}(count,4) = t(ii,jj);                    % theta of pixels
                    count = count + 1;
                end
            end
        end
        [~,Ipm] = max(power_peak1{i}(:,3));     % index of power-maximum
        wavelength_max(i,1) = power_peak1{i}(Ipm,2);
    end
    if peaks(i,2) ~= 0
        count = 1;
        for jj = 1:fz
            for ii = 1:fz
                if (t(ii,jj) >= peaks(i,2)-(step/2)) && (t(ii,jj) <= peaks(i,2)+(step/2)) ...
                            && (r(ii,jj) >= rmin) && (r(ii,jj) <= rmax)
                    power_peak2{i}(count,1) = r(ii,jj);
                    power_peak2{i}(count,2) = fz/power_peak2{i}(count,1);
                    power_peak2{i}(count,3) = abs(frameFFT2{i}(ii,jj));
                    power_peak2{i}(count,4) = t(ii,jj);
                    count = count + 1;
                end
            end
        end
        [~,Ipm2] = max(power_peak2{i}(:,3));
        wavelength_max(i,2) = power_peak2{i}(Ipm2,2);
    end
    if peaks(i,3) ~= 0
        count = 1;
        for jj = 1:fz 
            for ii = 1:fz
                if (t(ii,jj) >= peaks(i,3)-(step/2)) && (t(ii,jj) <= peaks(i,3)+(step/2)) ...
                            && (r(ii,jj) >= rmin) && (r(ii,jj) <= rmax)
                    power_peak3{i}(count,1) = r(ii,jj);
                    power_peak3{i}(count,2) = fz/power_peak3{i}(count,1);
                    power_peak3{i}(count,3) = abs(frameFFT2{i}(ii,jj));
                    power_peak3{i}(count,4) = t(ii,jj); 
                    count = count + 1;
                end
            end
        end
        [~,Ipm3] = max(power_peak3{i}(:,3));
        wavelength_max(i,3) = power_peak3{i}(Ipm3,2);
    end
    if peaks(i,4) ~= 0
        count = 1;
        for jj = 1:fz
            for ii = 1:fz
                if (t(ii,jj) >= peaks(i,4)-(step/2)) && (t(ii,jj) <= peaks(i,4)+(step/2)) ...
                            && (r(ii,jj) >= rmin) && (r(ii,jj) <= rmax)
                    power_peak4{i}(count,1) = r(ii,jj);
                    power_peak4{i}(count,2) = fz/power_peak4{i}(count,1);
                    power_peak4{i}(count,3) = abs(frameFFT2{i}(ii,jj));
                    power_peak4{i}(count,4) = t(ii,jj);
                    count = count + 1;
                end
            end
        end
        [~,Ipm4] = max(power_peak4{i}(:,3));
        wavelength_max(i,4) = power_peak4{i}(Ipm4,2);
    end
    if ~isempty(peak_ft{i}) 
        count = 1; 
        for jj = 1:fz
            for ii = 1:fz
                if (t(ii,jj) >= peak_ft{i}-(step/2)) && (t(ii,jj) <= peak_ft{i}+(step/2)) ...
                            && (r(ii,jj) >= rmin) && (r(ii,jj) <= rmax)
                    power_peak_ft{i}(count,1) = r(ii,jj);
                    power_peak_ft{i}(count,2) = fz/power_peak_ft{i}(count,1);
                    power_peak_ft{i}(count,3) = abs(frameFFT2{i}(ii,jj));
                    power_peak_ft{i}(count,4) = t(ii,jj);
                    count = count + 1;
                end
            end
        end
        [~,Ipm_ft] = max(power_peak_ft{i}(:,3));
        wavelength_ft(i,1) = power_peak_ft{i}(Ipm_ft,2);
    end
end
clear('Ipm','Ipm2','Ipm3','Ipm4','Ipm_ft')
    
%% dispfig == -1 (maps of backtransformations)

if fz == 101

    no_peak = zeros(51); no_peak(2:end-1,2:end-1) = 1;
                
    map_max = ones(size(image,1),size(image,2));
    f = 1;
    for v = 1:0.5:nframes_v
        for h = 1:0.5:nframes_h
            updown = ((v-1)*(fz-1))+25:((v-1)*(fz-1))+75;
            leftright = ((h-1)*(fz-1))+25:((h-1)*(fz-1))+75;
            if ~isequal(backtransf_max{f}(2,2),1)
                map_max(updown,leftright) = backtransf_max{f}(25:75,25:75);
            else
                map_max(updown,leftright) = no_peak;
            end
            f = f+1;
        end
    end

    map_ft = ones(size(image,1),size(image,2));
    f = 1;
    for v = 1:0.5:nframes_v
        for h = 1:0.5:nframes_h
            updown = ((v-1)*(fz-1))+25:((v-1)*(fz-1))+75;
            leftright = ((h-1)*(fz-1))+25:((h-1)*(fz-1))+75;
            if ~isequal(backtransf_ft{f}(2,2),1)
                map_ft(updown,leftright) = backtransf_ft{f}(25:75,25:75);
            else
                map_ft(updown,leftright) = no_peak;
            end
            f = f+1;
        end
    end
    
    if dispfig == -1
        figure; imshow(originalimage,[]); title('original image');
        figure; imshow(real(map_max),[-0.08 0.1]); % display range adjusted for visibility
        title('map of backtransformations along direction of max peak')
        figure; imshow(map_ft,[-0.08 0.1]); % display range adjusted for visibility
        title(['map of backtransformations along direction of features (\theta = ' ...
            num2str(tmin), ' to ' num2str(tmax), ' )'])
    end
    
elseif (dispfig == -1) && (fz ~= 101) && (fz ~= 51)
    disp('Sorry, this function currently only produces maps for a framesize of 51 or 101 px')
    
else
    map_max = ones(size(image,1),size(image,2));  %#ok<PREALL>
    map_ft = ones(size(image,1),size(image,2)); %#ok<PREALL>
    
end

if fz == 51
    nframes_h51 = round((size(image,2)-1)/(fz-1)) ; % horizontal frames, non-overlapping
    nframes_v51 = round((size(image,1)-1)/(fz-1)) ; % vertical frames, non-overlapping
    
    no_peak = zeros(27); no_peak(2:end-1,2:end-1) = 1;
                
    map_max = ones(size(image,1),size(image,2));
    f = 1;
    for v = 1:0.5:nframes_v51
        for h = 1:0.5:nframes_h51
            updown = ((v-1)*(fz-1))+13:((v-1)*(fz-1))+39;
            leftright = ((h-1)*(fz-1))+13:((h-1)*(fz-1))+39;
            if ~isequal(backtransf_max{f}(2,2),1)
                map_max(updown,leftright) = backtransf_max{f}(13:39,13:39);
            else
                map_max(updown,leftright) = no_peak;
            end
            f = f+1;
        end
    end

    map_ft = ones(size(image,1),size(image,2));
    f = 1;
    for v = 1:0.5:nframes_v51
        for h = 1:0.5:nframes_h51
            updown = ((v-1)*(fz-1))+13:((v-1)*(fz-1))+39;
            leftright = ((h-1)*(fz-1))+13:((h-1)*(fz-1))+39;
            if ~isequal(backtransf_ft{f}(2,2),1)
                map_ft(updown,leftright) = backtransf_ft{f}(13:39,13:39);
            else
                map_ft(updown,leftright) = no_peak;
            end
            f = f+1;
        end
    end
    
    if dispfig == -1
        figure; imshow(originalimage,[]); title('original image');
        figure; imshow(map_max,[-0.08 0.1]); % display range adjusted for visibility
        title('map of backtransformations along direction of max peak')
        figure; imshow(map_ft,[-0.08 0.1]); % display range adjusted for visibility
        title(['map of backtransformations along direction of features (\theta = ' ...
            num2str(tmin), ' to ' num2str(tmax), ' )'])
    end
    
else
    map_max = ones(size(image,1),size(image,2)); 
    map_ft = ones(size(image,1),size(image,2));
end

%% dispfig == 0.1 (text output)
peak_ft_nonempty = cell2mat(peak_ft(~cellfun('isempty',peak_ft)));   
    % (i.e. do not apply the following to empty cells)

if isequal(dispfig,0.1)
    fprintf(['\nThe image is of size ',num2str(size(image,1)),' x ',...
        num2str(size(image,2)),' px and was split into ',num2str(nframes),...
        ' frames of size ',num2str(fz),' x ',num2str(fz),' px.\n']);
    fprintf(['Parameters used: Sector width = ',num2str(step),...
        '°. Sensitive radii = ',num2str(rmin),'-',num2str(rmax),...
        ' px. Threshold = ',num2str(threshold),'x mean intensity.\n']);
    fprintf(['Acceptable angles (theta) for features: ',num2str(tmin),...
        '° to ',num2str(tmax),'°. Features detected between ',...
        num2str(min(peak_ft_nonempty)),'° and ',num2str(max(peak_ft_nonempty)),'°.\n']);
end

%% dispfig >= 1 (multiplot for specific frame number)

if dispfig >= 1
    x = dispfig;    
    figure; set(gcf,'color','w');
%1
    subplot(2,5,1); imshow(frames{x},[]); title(['frame number: ', num2str(x) ]); hold on
        text(0.12*fz,1.04*fz,['frame size: ' num2str(fz) ' x ' num2str(fz) ' px'],...
        'Color',[0 0 0],'Fontsize',8)
%2
    subplot(2,5,2); imshow(frames{x}.*window);title('frame after windowing'); hold on
        text(0.02*fz,1.04*fz,['Gaussian window mask {\sigma}=' num2str(0.4*fz,2) ' px'],...
        'Color',[0 0 0],'Fontsize',8)
%3  
    % FFT is displayed logarithmically for clarity
    subplot(2,5,3); imshow(log(1+abs(frameFFT2{x})),[]); 
    title('inspected area in the FFT2'); hold on;
        for e = 1:360
            plot(c + sind(e)*rmax, c + cosd(e)*rmax,'.','Markersize',(fz/rmax),...
                'Color',[1. 1. .99]); hold on
            plot(c + sind(e)*rmin, c + cosd(e)*rmin,'.','Markersize',(fz/rmax),...
                'Color',[1. 1. .99]); hold on
        end
        txt = {['sector width ' num2str(step) '° (50% overlap)'],...
                ['sensitive radius ' num2str(rmin) '-' num2str(rmax) ' px']};
        text([0.06 0.06]*fz,[1.04 1.10]*fz,txt,'Color',[0 0 0],'Fontsize',8); hold off
%4
    subplot(2,5,4); 
    imshow(backtransf_max{x},[]); title('signal of max. intensity'); hold on;
    if ~isequal(peaks(x,1),0)
        text(0.12*fz,1.04*fz,['modes with {\theta} = ' ...
            num2str(peaks(x,1)) '° ± ' num2str(step/2) '°'],'Color',[0 0 0],'Fontsize',8);
    else
        text(c-20,c,['no lineaments' newline 'detected in image'],...
            'Color',[0 0 0],'Fontsize',10);
    end
    hold off;

%5
    subplot(2,5,5); 
    imshow(backtransf_max_alt{x},[]); title('signal of max. intensity'); hold on;
    if ~isequal(peaks(x,1),0)
        txt = {['modes with {\theta} = ' num2str(peaks(x,1)) '° ± ' num2str(step/2) '°'],...
            'and rmin \leq r \leq rmax'};
        text([0.12 0.16]*fz,[1.04 1.11]*fz,txt,'Color',[0 0 0],'Fontsize',8);
    else
        text(c-20,c,['no lineaments' newline 'detected in image'],...
            'Color',[0 0 0],'Fontsize',10);
    end
%     imshow(backtransf_ft{x},[]); title('direction of features'); hold on;
%     if ~isempty(peak_ft{x})
%         text(0.22*fz,1.04*fz,['(for modes within [rmin rmax])' newline 'direction is {\theta}=' num2str(peak_ft{x}) '°'],...
%             'Color',[0 0 0],'Fontsize',8)
%     else
%         text(c-15,c,['no features' newline 'detected within' newline 'specified directions'],...
%             'Color',[0 0 0],'Fontsize',10);
%     end

% 6 (second row)
    subplot(2,2,3); ax1 = gca;
        stairs(sectors(1:nsec+1),sectorsum(x,1:nsec+1),'Color','k','parent',ax1); hold on; 
        avg = plot(xlim,[mean(sectorsum(x,:)) mean(sectorsum(x,:))],'-b'); hold on; 
        if ~isempty(peaks)
          peakmarks = plot(peaks(x,:),intensity(x,:)+max(sectorsum(x))/40,'v',...
              'MarkerSize',6,'Color',[0 .7 .7],'MarkerFaceColor',[0 .7 .7]); hold on;                                                                   
          text(peaks(x,:)+2.5,intensity(x,:)+(max(sectorsum(x))/40),...
              num2str((1:numel(intensity(x,:)))')); 
          hold off;
        legend([avg peakmarks],...
          {['mean' newline 'intensity of' newline 'all sectors'],...
          ['{}' newline 'peaks' newline '(sorted by' newline 'descending' newline ...
          'height)']},...
          'Location','northeastoutside');  
        else
            legend(avg,{['mean' newline 'intensity of' newline 'all sectors']},...
            'Location','northeastoutside'); 
        end
        title({['cumulative FFT-intensities per sector', newline '{}' newline '{}']});
        ax1 = gca; 
            set(ax1,'Position',[0.1300 0.1100 0.2685 0.3412],...
                'XTick',[0:30:180],'XTicklabel',[0:30:180]);%#ok<NBRAK>
            xlabel(ax1,'{\theta} : direction in the Fourier space [°]'); 
            ylabel(ax1,'\Sigma of FFT intensities')
            axis(ax1,[-1 183 0 inf]); %up to 28500
        ax2 = axes('Position',[0.1300 0.1100 0.2685 0.3412]);
            axis(ax2,[-91 93 0 inf]);
            set(ax2,'Color','none','XAxislocation','top','YAxislocation','right',...
                'XTick',[-90:30:90],'XTicklabel',[-90:30:90],...
                'YTick',[],'YTicklabel',[]); %#ok<NBRAK>
            xlabel(ax2,'{\alpha} : direction in the image space [°]');

% 7
    subplot(2,2,4); 
    if ~isempty(peak_ft{x})
        lay_id_text = ['(features = peak no.' num2str(peak_ft_index{x}) ')'];
    else
        lay_id_text = '(no features detected)';
    end

    if peaks(x,1) ~= 0
        plot(power_peak1{x}(:,2).*resolution,power_peak1{x}(:,3),'-','Color','k',...
            'Marker','d','MarkerFaceColor','w','Markersize',6); hold on;
        legend({['Peak 1: max. at ' num2str(wavelength_max(x,1),2) ...  
            ' px, (' num2str(wavelength_max(x,1).*resolution,2) ' m)' ...
            newline lay_id_text]},'Location','northwest');
    end

    if peaks(x,2) ~= 0
        plot(power_peak2{x}(:,2).*resolution,power_peak2{x}(:,3),'-','Color',[1 0 0],...
            'Marker','s','MarkerFaceColor',[1 0 0],'Markersize',6);hold on;
        legend({['Peak 1: max. at ' num2str(wavelength_max(x,1),2) ...
            ' px, (' num2str(wavelength_max(x,1).*resolution,2) ' m)'] ... 
                ['Peak 2: max. at ' num2str(wavelength_max(x,2),2) ...
            ' px, (' num2str(wavelength_max(x,2).*resolution,2) ' m)' ... 
            newline lay_id_text]},'Location','northwest');
    end

    if peaks(x,3) ~= 0
         plot(power_peak3{x}(:,2).*resolution,power_peak3{x}(:,3),'-','Color','b',...
             'Marker','o','MarkerFaceColor','w','Markersize',6); hold on;
         legend({['Peak 1: max. at ' num2str(wavelength_max(x,1),2) ...
            ' px, (' num2str(wavelength_max(x,1).*resolution,2) ' m)'] ... 
               ['Peak 2: max. at ' num2str(wavelength_max(x,2),2) ...
            ' px, (' num2str(wavelength_max(x,2).*resolution,2) ' m)'] ...
               ['Peak 3: max. at ' num2str(wavelength_max(x,3),2) ...
            ' px, (' num2str(wavelength_max(x,3).*resolution,2) ' m)' ...
            newline lay_id_text]},'Location','northwest');
    end

    if peaks(x,4) ~= 0
        plot(power_peak4{x}(:,2).*resolution,power_peak4{x}(:,3),'-','Color','m',...
            'Marker','^','MarkerFaceColor','m','Markersize',6); hold on;
        legend({['Peak 1: max. at ' num2str(wavelength_max(x,1),2) ...
            ' px, (' num2str(wavelength_max(x,1).*resolution,2) ' m)'] ... 
                ['Peak 2: max. at ' num2str(wavelength_max(x,2),2) ...
            ' px, (' num2str(wavelength_max(x,2).*resolution,2) ' m)'] ...
                ['Peak 3: max. at ' num2str(wavelength_max(x,3),2) ...
            ' px, (' num2str(wavelength_max(x,3).*resolution,2) ' m)'] ...
                ['Peak 4: max. at ' num2str(wavelength_max(x,4),2) ...
            ' px, (' num2str(wavelength_max(x,4).*resolution,2) ' m)' ...
            newline lay_id_text]},'Location','northwest');
    end
    title({['FFT power spectrum along the direction of intensity peaks', newline '{}']});
    xlabel('wavelength of the Fourier mode [in m]');
    ylabel('intensity of the Fourier mode');
    yticklabels({0,1,2,3,4,5})
    xlim([4 inf])

end

end % end of function       
                
