clear
load 128_90_250_SystemMatrix.mat
imageResolution = 128;
%imageResolution = 64;
original_image = phantom(imageResolution);
A = SystemMatrix;
w = 1;
x=original_image';
x=x(:);
y = full(SystemMatrix) * x;
x0 = zeros(128^2,1);
t = 100;
x = sart( A, w, y, x0, t);
reconstructed_image = reshape(x,128,128);
reconstructed_image=reconstructed_image';
imshow(reconstructed_image)