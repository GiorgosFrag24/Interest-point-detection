clear all;
close all;
I = imread('sunflowers.png');
I = im2double(I);
I1 = rgb2gray(I);
info = size(I1);

sigma = 2;
p = 2.5;
s = 1.5;
N = 6;

interest_points1 = exer02_1(I1,info,sigma,p);
interest_points_visualization(I,interest_points1);
print -djpeg exer021.jpg

interest_points2 = exer02_2(I1,sigma,p,s,N,info);
interest_points_visualization(I,interest_points2);
print -djpeg exer022.jpg

interest_points3 = exer02_3(I1,sigma);
interest_points_visualization(I,interest_points3);
print -djpeg exer023.jpg

interest_points4 = exer02_4(I1,sigma,s,N);
interest_points_visualization(I,interest_points4);
print -djpeg exer024.jpg

interest_points5 = exer02_5(I1,sigma);
interest_points_visualization(I,interest_points5);
print -djpeg exer025.jpg

interest_points5b = exer02_5b(I1,sigma,s,N);
interest_points_visualization(I,interest_points5b);
print -djpeg exer025b.jpg