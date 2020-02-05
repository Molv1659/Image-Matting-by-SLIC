clear all
A = imread('.\data\3.jpg');
[L,N] = superpixels(A,500); 
% [L,N] = myslic(A,100000,5,7);

figure
BW = boundarymask(L);
imshow(BW);
figure
imshow(imoverlay(A,BW,'cyan'),'InitialMagnification',67)

%»­³ÉË®Ä«·ç
% outputImage = zeros(size(A),'like',A);
% idx = label2idx(L);
% numRows = size(A,1);
% numCols = size(A,2);
% for labelVal = 1:N
%     redIdx = idx{labelVal};
%     greenIdx = idx{labelVal}+numRows*numCols;
%     blueIdx = idx{labelVal}+2*numRows*numCols;
%     outputImage(redIdx) = mean(A(redIdx));
%     outputImage(greenIdx) = mean(A(greenIdx));
%     outputImage(blueIdx) = mean(A(blueIdx));
% end
% 
% figure
% imshow(outputImage,'InitialMagnification',67)
% imwrite(outputImage,'3.jpg')