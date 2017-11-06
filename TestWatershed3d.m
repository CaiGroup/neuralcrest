% the structure of the file name:
addpath(genpath('C:\Users\alignell\Desktop\20160411_UCSB\antti'));
nameDummy = 'RGB';
% ilastik stores results in a h5 container. Exported_data is like a folder
% containing the probability map for your image. Each pixel has number
% between 0 1
file = h5read([nameDummy,'_bin2_Probabilities.h5'],'/exported_data/');
% file is an array with NSpatialDim+NChannel+1 dimensions. the last channel contains the number 
% of classes used during training. (2 is default, 2nd considered
% foreground)
pred = squeeze(file(2,:,:,:));
%pred = squeeze(file(1,:,:,:));
% undo internal axis swapping of ilastik
pred = permute(pred,[2,1,3]);
%% Visualize summed prediction for foreground (Nuclei).
imshow(sum(pred,3),[])

%%  read the image; make sure membrane is first (second) channel
clear im;
p = 2;% Ratio between binning of ilastik segmentation image and actual image used for watershed.
for k = 1 : size(pred,3)%*p
    %temp = imread([nameDummy,'_bin6.tif'],k);
    temp = imread([nameDummy,'.tif'],k);
    im(:,:,k) = temp(:,:,1);
end

im = im(1:size(pred,1)*p,1:size(pred,2)*p,:);
%% blurred version of the membrane channel
im_blurr = imdilate(im,strel3D('sphere',3));

imshow(max(im_blurr,[],3),[])

%% Take ilastik prediction map and threshold
seg = pred>.995; % number between 0,1 (.6 one of the best values)
seg = bwareaopen(seg,150); % require minimum nucleus size
seg = seg - bwareaopen(seg,30000); % require maximum nucleus size
%% interpolate the image onto the bin2 image;
%seg = imresize(seg,3);
if p > 1
    [X,Y,Z] = meshgrid(1/p:1/p:size(pred,2),1/p:1/p:size(pred,1),1:size(pred,3));
    %[X,Y,Z] = meshgrid(1/p:1/p:size(pred,2),1/p:1/p:size(pred,1),1/p:1/p:size(pred,3));
    %[X1,Y1,Z1] = meshgrid(1:size(pred,2),1:size(pred,1),1:size(pred,3));
    seg2 = interp3(seg,X,Y,Z);
else
    seg2 = seg;
end

%% at all positions where seg2 == 1, imimposemin imposes a minimum on im_blurr, 
%, which we use as seed for watershed. 
I2 = imimposemin(im_blurr,seg2);

%% 
% use this iamge I2, with the seeds for watershed, and run watershed in 3D.
% This gives us a label matrix L
Label = watershed(I2);
%% get rid of the boundary touching cells;
Label2 = imclearborder(Label);

%% save label to disk; 
for k = 1 : size(Label,3)
     imwrite(Label(:,:,k),'Label_WatershedMatrix.tif','compression','none','WriteMode','append');
end
save('Label_WatershedMatrix.mat','Label');

%% try to find cell number x
ll= zeros(size(Label),'uint16');
for x = [183,272]
    
ll = ll+ uint16(Label == x);
imshow(max(ll,[],3),[])
end
%% Visually check the result of watershed segmentation slice by slice
delete('seg_RGB.tif')
for k = 1 : size(im,3)
    
    temp = zeros(size(im,1),size(im,2),3,'uint8');
    per = Label(:,:,k) == 0; % use LABEL to have all cells, including boundary
   % per = Label2(:,:,k) == 0; % only for perfeclty intact ones;
    
    temp(:,:,1) = im(:,:,k); % this is the memrbane channel
    temp(:,:,2) = uint8(per)*100; % this is the perimter obtained by LABEL or LABEL2
    %imshow(temp,[])
    %M(k) = getframe;
    % write an image to disc, that shows membrane in red, outline in green
    imwrite(temp,'seg_RGB.tif','compression','none','WriteMode','append');
end

%%  generate a colourmap, random. 
cmp = jet(double(max(Label(:))));
cmp = cmp(randperm(size(cmp,1)),:);
%% % generate a RGB Stack showing pseudo colour of each cell; Look at it in IMARIS
RGBFName = 'seg_label_RGB.tif';
delete(RGBFName);
R = 0*Label;
G = R;
B = R;
s = regionprops(Label,'Area');
volume = cat(1,s.Area);
vol_min = 500;
vol_max = 5000; % maximum volume of cell in pixels; may need to be changed
out = union(find(volume<vol_min),find(volume>vol_max));
cmp(out,:) = 0*cmp(out,:);
for k = 1 : size(Label,3)
    temp = label2rgb(Label(:,:,k),cmp,'k');
    R(:,:,k) = temp(:,:,1);
    G(:,:,k) = temp(:,:,2);
    B(:,:,k) = temp(:,:,3);
    
    imwrite(temp,RGBFName,'compression','none','WriteMode','append');
end
    
%% 
%
s = regionprops(Label);
area = cat(1,s.Area);
index = intersect(find(area>vol_min),find(area<vol_max));
colordef black;
hold on 
%% do the matlab way of showing data in 3D
for ki = 1:length(index);
    k = index(ki);
p = patch(isosurface(Label == k,.5));
 isonormals(Label,p)
 hold on;
 set(p,'FaceColor',cmp(k,:),'EdgeColor','none');
end