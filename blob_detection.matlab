function [interest_points3] = exer02_3(I,sigma)

%% Create the smoothed image
close all;
n = ceil(3*sigma)*2+1; 
Gs = fspecial('gaussian',n,sigma);                                          %create the Gs 2D smoothing kernel
Is = imfilter(I,Gs,'symmetric');                                            %equivalent to I*Gs 

%% 2.3.1 Compute the Hessian matrix (2nd degree grad) and accordingly the criterion R
[Lx, Ly] = gradient(Is);                                                    %compute 1st degree grads
[Lxx, Lxy] = gradient(Lx);                                                  %compute 2nd degree grads
[~, Lyy] = gradient(Ly);                                                    %compute 2nd degree grads

R = Lxx.*Lyy - Lxy.^2;

%% 2.3.2 Find points of interest
thcorn = 0.05;
Rmax = max(max(R));

ns = ceil(3*sigma)*2+1;
B_sq = strel('disk',ns);
Cond1 = ( R==imdilate(R,B_sq) );
[inti,intj] = find(Cond1 & R>thcorn*Rmax);
[x,~] = size(inti);
temp = ones(x,1)*sigma;

interest_points3 = [intj inti temp];
end

