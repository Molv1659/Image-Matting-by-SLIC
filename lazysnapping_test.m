clear all
RGB = imread('.\data\3.jpg');
figure; 
imshow(RGB)
h1 = impoly(gca,'Closed',false);
foresub = getPosition(h1);
foresub = int32(foresub);
foregroundInd = sub2ind(size(RGB),foresub(:,2),foresub(:,1));
figure; 
imshow(RGB)
h2 = impoly(gca,'Closed',false);
backsub = getPosition(h2);
backsub = int32(backsub);
backgroundInd = sub2ind(size(RGB),backsub(:,2),backsub(:,1));
[L,N] = myslic(RGB,300,5,10);
%L = superpixels(RGB,1500);
% BW = boundarymask(L);
% imshow(imoverlay(RGB,BW,'cyan'),'InitialMagnification',67)
BW = lazysnapping(RGB,L,foregroundInd,backgroundInd);
figure;
imshow(BW);
maskedImage = RGB;
maskedImage(repmat(~BW,[1 1 3])) = 0;
figure; 
imshow(maskedImage)


