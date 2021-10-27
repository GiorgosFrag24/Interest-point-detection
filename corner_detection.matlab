function [interest_points] = exer02_1(I,info,sigma,p)

%% 2.1.1 Compute stractural tensor J
close all;
n = ceil(3*sigma)*2+1; 
Gs = fspecial('gaussian',n,sigma);                                          %create the Gs 2D smoothing kernel
n = ceil(3*p)*2+1; 
Gp = fspecial('gaussian',n,p);                                              %create the Gp 2D smoothing kernel
Is = imfilter(I,Gs,'replicate');                                            %equivalent to I*Gs     

[Isx,Isy] = gradient(Is);                                           
J1 = imfilter((Isx.*Isx),Gp,'replicate');                                   %J1, equivalent to the given type
J2 = imfilter((Isx.*Isy),Gp,'replicate');                                   %J2, equivalent to the given type
J3 = imfilter((Isy.*Isy),Gp,'replicate');                                   %J3, equivalent to the given type

%% 2.1.2 Compute the eigenvalues L+, L- for the tensor J
Lplus = 0.5*(J1+J3+sqrt((J1-J3).^2+4*J2.^2));
Lminus = 0.5*(J1+J3-sqrt((J1-J3).^2+4*J2.^2));

%% 2.1.3 Compute the cornerness criterion, and accordingly the points of interest in the image
k = 0.05;
thcorn = 0.005;

R = Lminus.*Lplus - k*(Lminus+Lplus).^2;
Rmax = max(max(R));

ns = ceil(3*sigma)*2+1;
B_sq = strel('disk',ns);
Cond1 = ( R==imdilate(R,B_sq));
[inti,intj] = find(Cond1 & R>thcorn*Rmax);
[x,~] = size(inti);
temp = ones(x,1)*sigma;

interest_points = [intj inti temp];
end