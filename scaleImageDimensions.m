% We assume there are two types of images produced, with different
% resolution + orientation etc.
%AimedMag  = 0.227; % pixel size in microns used to analyse the data;
addpath(genpath('C:\Users\alignell\Desktop\20160411_UCSB\antti'))
AimedMag  = 0.1135; % pixel size in microns used to analyse the data;

nChan     = 6; % number of channels used
dapiChan  = 6; % specifies a nucleus labeled channel
hybnum = 1;
sample = '3b2_lower';

% scale images to the correct dimensions
scope = 'olympus';%'apotome';
imageName = '0025_pos3b2_lower_Ltip_4som_Cy7Foxd3_A647Sox9_A594Bcl2_Cy3bSip1_A488ActB_100x_MMStack_Pos0.ome.tif';

nonUniformIlluminationName = 'Illucorr_cy7_647_594_cy3b_488_Dapi_no2_MMStack_Pos0.ome.tif';

% apotopmeImage serves as the reference; Segmented image, that we produced
% watershed LABEL mask from; Always produce LABEL first. No Label, script
% not running.
nucleusImage = 'RGB.tif';
nucleusChan = 2; % the nucleus channel in the nucleus image; %%now using membrane channel for this due to better results

% settings of the two different types of images. These values might need to
% be changed, depending on the acquisition used.
if strcmp(scope,'olympus');
    mag = 100;
    pixelSize = 13; % in microns;
    dx = pixelSize/mag;
    dz = .5; % in microns
elseif strcmp(scope,'apotome');
    mag = 40;
    pixelSize = 4.54; % in microns
    dx = pixelSize/mag;
    dz = .5; % in microns
end
% determine the scaling factor.
scaleFactorX = 1; %AimedMag/dx;
%scaleFactorZ = AimedMag/dz;
scaleFactorZ = 1;

% read the data; This is then stored as a cell array with 1*nChan entries.
fileSize = length(imfinfo(imageName));
stackSize = fileSize/nChan;
im = cell(1,nChan);
im_bg = cell(1,nChan); % background cell array, for correction of non-uniform illumination;
for c = 1:nChan
    im_bg{c} = imread(nonUniformIlluminationName,c);
    for k = 1:stackSize       
        if c < 6
            im{c}(:,:,k) = (double(imread(imageName,c+(k-1)*nChan)));%./((double(im_bg{c}+1)/max(abs(double(im_bg{c}(:))))));
            %im{c}(:,:,k) = (1000*double(imread(imageName,c+(k-1)*nChan)))./double((im_bg{c}+1));
        else
            im{c}(:,:,k) = double(imread(imageName,c+(k-1)*nChan));
        end
    end
end
%% scale the image to the write Aimed PixelSize;
im_scaled = cell(size(im));
[X,Y] = meshgrid(1:scaleFactorX:size(im{1},1),1:scaleFactorX:size(im{1},2));
newNSlices = round(size(im{1},3)/scaleFactorZ);
for c = 1 : nChan
    im_scaled{c} = imresize(im{c},1/scaleFactorX);
    si = size(im_scaled{c});
    scaled = permute(im_scaled{c},[3,1,2]);
    scaled_scaled = imresize(scaled,[newNSlices,size(scaled,2)]);
    im_scaled{c} = ipermute(scaled_scaled,[3,1,2]);
end
%% write the scaled images to disc
for c = 1 : nChan-1
    
    for k = 1 : size(im_scaled{c},3)
        
        imwrite(uint16(im_scaled{c}(:,:,k)),sprintf('hyb%d_Pos%s_stack_chan_%d_illucorr.tif',hybnum,sample,c),'tiff','compression','none','WriteMode','append');
    end
end
%% Use this only if you want to use post modified im_scaled
im_scaled_process = cell(1,nChan);
for c = 1:nChan
    for k = 1:stackSize       
        if c < 6
            im_scaled_process{c}(:,:,k) = double(imread(sprintf('hyb%d_Pos%s_stack_chan_%d_illucorr.tif',hybnum,sample,c),k));
            %im{c}(:,:,k) = (1000*double(imread(imageName,c+(k-1)*nChan)))./double((im_bg{c}+1));
        else
            im_scaled_process{c}(:,:,k) = im_scaled{c}(:,:,k);
        end
    end
end
im_scaled = im_scaled_process;
%% remove background in dapi channel;
% c = 6;
% sD = strel3D('sphere',10);
% %dapi_bg = imerode(im_scaled{c},sD);
% %dapi_bg = imdilate(dapi_bg,sD);
% %im_scaled{c} = im_scaled{c}-dapi_bg;
% for k = 1 : size(im_scaled{c},3);
%     
%     tmp = im_scaled{c}(:,:,k);
%     tmp = imopen(tmp,strel('disk',60));
%     im_scaled{c}(:,:,k) = im_scaled{c}(:,:,k)-tmp;
% end

%%  Read the image with the nucleus and the membrane (The one used for watershed seg.)
NucleusStackSize = length(imfinfo(nucleusImage));

imNuc = cell(1,2);

for k = 1 : NucleusStackSize
    temp = imread(nucleusImage,k);
    imNuc{1}(:,:,k) = temp(:,:,1);
    imNuc{2}(:,:,k) = temp(:,:,2);
end

%% manual rotation from im_scaled into imNuc frame of reference; 
% determine the rotation angle; 
chan = 5;
angle = 0;
im_rot = imrotate(im_scaled{chan},angle,'crop');
%imshow(im_rot(:,:,30),[])

% Determine drift between images; 
% beware that the hotspots images are smaller than the nucleus iamge; so we
% crop the nucleus image; 
%xcrop = 1000:2172;
%ycrop = 130:1302;
% correlate nucleus image with one of the channels to determine shift
% between scopes
%image1 = imdilate(imNuc{nucleusChan}(ycrop,xcrop,:),strel3D('sphere',3));
%image2 = imdilate(uint8(20sixe(0*mat2gray(im_scaled{6})),strel3D('sphere',3));
%[shift,image1,image2] = xcorr3fft(image1,image2);
%[shift,image1,image2] = xcorr3fft(imNuc{nucleusChan}(ycrop,xcrop,:),uint8(255*mat2gray(im_rot)));
[shift,image1,image2] = xcorr3fft(imNuc{nucleusChan}(:,:,:),uint8(255*mat2gray(im_rot)));
% retruns image1 = imNuc{nucleusChan}(ycrop,xcrop,:)
% image2 = uint8(255*mat2gray(im_scaled{3}))
% shift = Translation from iamge2 into image1, a vector with 3 components. 
%%
 % takes image2, and shifts it according to shift. 
close all; 
%figure, imshow(max(image1,[],3),[])
image2_shift = 0*image2;
if shift(1)<0 && shift(2)>0
    image2_shift(-shift(1)+(1:(size(image2,1)+shift(1))),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
        image2(1:(size(image2,1)+shift(1)),1:size(image2,2)-shift(2),:);
elseif shift(1)>0 && shift(2)<0
    image2_shift((1:size(image2,1)-shift(1)),1:size(image2,2)+shift(2),:) = ...
        image2((shift(1)+1):size(image2,1),-shift(2)+(1:size(image2,2)+shift(2)),:);
elseif shift(1)<=0 && shift(2)<=0
    image2_shift(-shift(1)+(1:size(image2,1)+shift(1)),1:size(image2,2)+shift(2),:) = ...
        image2((1:size(image2,1)+shift(1)),-shift(2)+(1:size(image2,2)+shift(2)),:);
elseif shift(1)>=0 && shift(2)>=0
    image2_shift((1:size(image2,1)-shift(1)),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
         image2((shift(1)+1):size(image2,1),1:size(image2,2)-shift(2),:);
elseif shift(1) == 0 && shift(2) == 0
    image2_shift = image2;
end


% visualize the result of the shift; 
% Red channel shows cropped field of view from apotome
% Green channel shows shifted version of image2 (containing hotspots);
%figure, imshow(max(image2_shift,[],3),[]);
tmp = zeros(size(image1,1),size(image1,2),3,'double');
tmp(:,:,1) = double(sum(image1,3))/size(image1,3); % mean intensity projection. 
tmp(:,:,2) = (sum(double(image2_shift)/size(image2,3),3))*8; % 5 is a scaling factor. 
tmp = uint8(255*mat2gray(tmp));
%tmp = imadjust(tmp,[.2 .3 0; .6 .7 1],[]);
imshow(tmp);

%% now inspect in a single slice cross section
k = 350;
figure,
tmp = zeros(size(image1,2),size(image1,3),3,'double');
tmp(:,:,1) = squeeze(image1(k,:,:));
tmp(:,:,2) = squeeze(image2_shift(k,:,:))*5;
tmp = uint8(255*mat2gray(tmp));
%tmp = uint8(tmp);
imshow(tmp);

% save('0026_values.mat','shift');



%%
%image1 = imNuc{nucleusChan}(ycrop,xcrop,:);
image1 = imNuc{1}(:,:,:);

%% summarize the dots in each cell, found using waterhed label. 
% First segment the hotspots using ilastk.
% after segmenting with ilastik, we have h5 files, that can be used to count the dots. 
% Count the dots in each cell;
points = struct('volume',[],'centroid',[],'MeanIntensity',[],'Intensity',[],...
    'NSpots',[],'NSpotsinVol',[],'NSpots_Cleared',[],'NSpotsinVol_Cleared',[]);
%% 
% summarizes information about cells, such as Volume, centre of mass
% position (centroid), Mean Intensity in the segmentation of the hotspots,
% total intensity; To get the total number of spots multiply by volume;
for chan = 1%:nChan-1
file = h5read(sprintf('hyb%d_Pos%s_stack_chan_%d_illucorr_Probabilities.h5',hybnum,sample,chan),'/exported_data/');
pred = squeeze(file(2,:,:,:)); % second class in ilastik is what we want. 
pred = permute(pred,[2,1,3]);
pt = .5;
image2 = pred>pt; % threshold value. can be low, typical value .5 or .3;

% rotate image 2, just as we rotated the reference channel above; 
image2 = imrotate(image2,angle,'crop'); % rotates the segmentation;
% shift the segmentation;
si1 = size(image1);
si2 = size(image1);
desiredSize = max(si1,si2);
padded_image = zeros(desiredSize,'uint8');
padded_image(1:size(image2,1),1:size(image2,2),1:size(image2,3)) = image2;
image2 = padded_image;
image2_shift = 0*image2;
if shift(1)<0 && shift(2)>0
    image2_shift(-shift(1)+(1:(size(image2,1)+shift(1))),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
        image2(1:(size(image2,1)+shift(1)),1:size(image2,2)-shift(2),:);
elseif shift(1)>0 && shift(2)<0
    image2_shift((1:size(image2,1)-shift(1)),1:size(image2,2)+shift(2),:) = ...
        image2((shift(1)+1):size(image2,1),-shift(2)+(1:size(image2,2)+shift(2)),:);
elseif shift(1)<=0 && shift(2)<=0
    image2_shift(-shift(1)+(1:size(image2,1)+shift(1)),1:size(image2,2)+shift(2),:) = ...
        image2((1:size(image2,1)+shift(1)),-shift(2)+(1:size(image2,2)+shift(2)),:);
elseif shift(1)>=0 && shift(2)>=0
    image2_shift((1:size(image2,1)-shift(1)),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
         image2((shift(1)+1):size(image2,1),1:size(image2,2)-shift(2),:);
elseif shift(1) == 0 && shift(2) == 0
    image2_shift = image2;
end


tmp = zeros(size(image1,1),size(image1,2),3,'double');
tmp(:,:,1) = double(sum(image1,3))/size(image1,3);
tmp(:,:,2) = (sum(double(mat2gray(image2_shift)*255)/size(image2,3),3))*4;
tmp = uint8(tmp);
figure, 
imshow(tmp);

% Label is determined in the script TestWatershed3D.m -> May make this nicer.
%s = regionprops(Label(ycrop,xcrop,:),image2_shift,'centroid','area','MeanIntensity','PixelValues');
s = regionprops(Label(:,:,:),image2_shift,'centroid','area','MeanIntensity','PixelValues');

points(chan).volume = cat(1,s.Area);
points(chan).centroid = cat(1,s.Centroid);

points(chan).MeanIntensity = cat(1,s.MeanIntensity);
points(chan).Intensity = points(chan).MeanIntensity.*points(chan).volume;


%s2 = regionprops(Label(ycrop,xcrop,:),im_scaled{chan},'MeanIntensity');
s2 = regionprops(Label(:,:,:),im_scaled{chan},'MeanIntensity');


points(chan).MeanDotIntensity = cat(1,s2.MeanIntensity);

bw = bwlabeln(image2_shift); % in the thresholded ilastik prediction map, which was rotated and shifted,
% we now look for connected groups of pixels. Each connected group gets a
% unique label, which we interpret as a hotspot. 
%s = regionprops(Label(ycrop,xcrop,:),bw,'centroid','area','MeanIntensity','PixelValues');
s = regionprops(Label(:,:,:),bw,'centroid','area','MeanIntensity','PixelValues');
NSpots = zeros(size(s));
for k = 1 : length(s);
    temp = s(k).PixelValues; % what labels are around in the cell k?
    temp = unique(temp); % extract the unique label values
    NSpots(k) = length(temp);%-1; % remove contribution from 0. 
end
points(chan).NSpots  = NSpots;
points(chan).NSpotsinVol  = NSpots./cat(1,s.Area);

%tmp = max(im_scaled{chan},[],3); % goal: Estimate the biggest background, and subtract!
%bg = imopen(tmp,strel('disk',10));
tmp = im_scaled{chan};

% rotate the originoal image 
tmp = imrotate(tmp,angle,'crop');
% shift hte original image; 
% shift the segmentation;
si1 = size(image1);
si2 = size(image1);
desiredSize = max(si1,si2);
padded_image = zeros(desiredSize,'uint16');
padded_image(1:size(image2,1),1:size(image2,2),1:size(image2,3)) = tmp;
tmp = padded_image;
tmp_shift = 0*tmp;
if shift(1)<0 && shift(2)>0
    tmp_shift(-shift(1)+(1:(size(image2,1)+shift(1))),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
        tmp(1:(size(image2,1)+shift(1)),1:size(image2,2)-shift(2),:);
elseif shift(1)>0 && shift(2)<0
    tmp_shift((1:size(image2,1)-shift(1)),1:size(image2,2)+shift(2),:) = ...
        tmp((shift(1)+1):size(image2,1),-shift(2)+(1:size(image2,2)+shift(2)),:);
elseif shift(1)<=0 && shift(2)<=0
    tmp_shift(-shift(1)+(1:size(image2,1)+shift(1)),1:size(image2,2)+shift(2),:) = ...
        tmp((1:size(image2,1)+shift(1)),-shift(2)+(1:size(image2,2)+shift(2)),:);
elseif shift(1)>=0 && shift(2)>=0
    tmp_shift((1:size(image2,1)-shift(1)),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
         tmp((shift(1)+1):size(image2,1),1:size(image2,2)-shift(2),:);
elseif shift(1) == 0 && shift(2) == 0
    tmp_shift = image2;
end



% for k = 1 : size(tmp,3)
%     bg = imopen(tmp(:,:,k),strel('disk',10));
%    tmp(:,:,k) = tmp(:,:,k)-bg; 
% end
bw = bwlabeln(image2_shift,6);
s2 = regionprops(bw,tmp_shift,'area','MeanIntensity','centroid');
DotVolumes = cat(1,s2.Area);
DotCentre = cat(1,s2.Centroid);
MeanDotIntensity = cat(1,s2.MeanIntensity);
% normalize MeanDotIntensityDistribution; 
[c,x] = ecdf(MeanDotIntensity);
%minimum = 1;%find(c<.1,1,'first');
maximum = find(c>.99,1,'first'); % normalize to 99 percent of the distribution
mI = MeanDotIntensity-min(MeanDotIntensity)+eps; % subtract minimum 
mI = mI/x(maximum)-eps/2; % divide by maximum % should give distribtion from 0 : 1+eps
points(chan).DotVolumes = DotVolumes;
points(chan).MeanIntensityOfDots = MeanDotIntensity;
points(chan).MeanIntensityOfDotsNormalized = mI;
figure, 
hist(mI,500);
ind = intersect(find(mI>.1),find(mI<1.5));% range, cutoff percent brightest values;
points(chan).index_filter = ind;
% display the selected group as red. 
figure,imshow(max(tmp_shift,[],3),[])
%figure, imshow(sum(tmp_shift,3),[])
hold on 
plot(DotCentre(:,1),DotCentre(:,2),'.') % All dots shown in blue
plot(DotCentre(ind,1),DotCentre(ind,2),'r.') % Selected dots shown in red. 

% now use the ind to determine which of the hot spots are in which cells; 
% to achieve this, we will remove all the labels of the bw matrix, that are
% not member of ind ; Make it a label matrix again, call bw2;
bw2 = bwlabeln(ismember(bw, ind),6);


 % in the thresholded ilastik prediction map, which was rotated and shifted,
% we now look for connected groups of pixels. Each connected group gets a
% unique label, which we interpret as a hotspot. 
%s = regionprops(Label(ycrop,xcrop,:),bw2,'centroid','area','MeanIntensity','PixelValues');
s = regionprops(Label(:,:,:),bw2,'centroid','area','MeanIntensity','PixelValues');
NSpots = zeros(size(s));
for k = 1 : length(s);
    temp = s(k).PixelValues; % what labels are around in the cell k?
    temp = unique(temp); % extract the unique label values
    NSpots(k) = length(temp);%-1; % remove contribution from 0. 
end
points(chan).NSpots_Cleared  = NSpots;
points(chan).NSpotsinVol_Cleared  = NSpots./cat(1,s.Area);

end
%pointssave = sprintf('hyb%d_points.mat',hybnum);
%save(pointssave,'points');
%% make some pretty images of the label matrix L;
% visualize the region of interst as 3D matrix, for inspection with e.g.
% Amira
bw = Label(ycrop,xcrop,:);
bw(bw == 1)  = 0;
cmp = hsv(double(max(bw(:))));
cmp = cmp(randperm(size(cmp,1)),:);

for k = 1 : size(bw,3)
    
    temp = label2rgb(bw(:,:,k),cmp,'k');
    %imshow(temp);
    imwrite(temp,'cc3d.tif','tiff','compression','none','WriteMode','append');
end

%% double check that MeanDotIntensity is reasonable.
figure, 
chan = 3;

mI1 = points(chan).MeanDotIntensity;
mI2 = points(chan).NSpotsinVol;

mI = mI2; % will display the mean intensity 

% find upper bound and lower bound; 
[c,x] = ecdf(mI);
minimum = 1;%find(c<.1,1,'first');
maximum = find(c>.95,1,'first');
mI = mI-min(mI)+eps;
mI = mI/x(maximum)-eps/2;
edges = 0:.1:1.5;
cmp = jet(length(edges));
[~,bin] = histc(mI,edges);
imshow(max(im_scaled{chan},[],3),[])
hold on
for i = 1 : length(edges)
    
    ind = find(bin == i-1);
    plot(points(chan).centroid(ind,1),points(chan).centroid(ind,2),'o',...
        'color',cmp(i,:),'LineWidth',2);
end
%% 

% find upper bound and lower bound; 
[c,x] = ecdf(mI1);
minimum = 1;%find(c<.1,1,'first');
maximum = find(c>.95,1,'first');
mI1 = mI1-min(mI1)+eps;
mI1 = mI1/x(maximum)-eps/2;



% find upper bound and lower bound; 
[c,x] = ecdf(mI2);
minimum = 1;%find(c<.1,1,'first');
maximum = find(c>.95,1,'first');
mI2 = mI2-min(mI2)+eps;
mI2 = mI2/x(maximum)-eps/2;

plot(mI1,mI2,'o'), axis([0 1.5 0 1.5])
xlabel('MeanDotIntensity normalized')
ylabel('SpotsinVolume normalized')
%% 
figure, 
diff = abs(mI1-mI2);
[~,bin] = histc(diff,edges);
imshow(max(im_scaled{chan},[],3),[])
hold on
for i = 1 : length(edges)
    
    ind = find(bin == i-1);
    plot(points(chan).centroid(ind,1),points(chan).centroid(ind,2),'o',...
        'color',cmp(i,:),'LineWidth',2);
end
