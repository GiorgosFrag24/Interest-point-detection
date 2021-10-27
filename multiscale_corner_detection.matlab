function [interest_points] = exer02_2(I,sigma0,p0,s,N,info)

%% 2.2.1 Apply the algorithm for different scales
for i=1:N
    sigmai(i) = (s^(i-1))*sigma0;                           %compute sigmai and pi for each scale
    pi(i) = (s^(i-1))*p0;
    temp(:,:,i) = {exer02_1(I,info,sigmai(i),pi(i))};       %and store the results for the algorithm in cell format
end

count = 1;                                                  %counter to store all interest points in one matrix
for i=1:N
    [sizex,~] = size(temp(:,:,i));
    for j=1:sizex
        points(count,:) = temp(j,:,i);                      %loop to do the transformation between cell format
        count = count+1;                                    %and the interest points matrix
    end
end

points = cell2mat(points);                                  %transformation complete

%% 2.2.2 Pick the appropriate scale for each point of interest
for i=1:N
    n = ceil(3*sigmai(i))*2+1; 
    Gs = fspecial('gaussian',n,sigmai(i));                  %create the Gs 2D smoothing kernel
    Is = imfilter(I,Gs,'symmetric');                        %equivalent to I*Gs 

    [Lx, Ly] = gradient(Is);                                %compute 1st degree grads
    [Lxx, ~] = gradient(Lx);                                %compute 2nd degree grads
    [~, Lyy] = gradient(Ly);                                %compute 2nd degree grads
    
    LOG(:,:,i) = {sigmai(i)^2.*abs(Lxx+Lyy)};               %compute LOG for each scale
end

LOG = cell2mat(LOG);                                        %convert cell type to matrix for use

[sizex,~] = size(points);                                   %how many points of interest we have to check?
count = 1;
for i=1:sizex
    [~,d] = find(sigmai == points(i,3));                    %find the scale to which the point under examination corresponds                 
    if (d==1)                                               %if its the first then check only with next
        if (LOG(points(i,2),points(i,1),1)>LOG(points(i,2),points(i,1),2))  
            interest_points(count,:) = points(i,:);         %and keep the point if it results in bigger metric LOG than it does for the next scale
            count = count+1;
        end
    elseif (d==N)                                           %if it is the last then check only with previous
       if (LOG(points(i,2),points(i,1),N-1)>LOG(points(i,2),points(i,1),N))
            interest_points(count,:) = points(i,:);
            count = count+1;
       end
    else
        if ((LOG(points(i,2),points(i,1),d)>LOG(points(i,2),points(i,1),d-1))&&(LOG(points(i,2),points(i,1),d)>LOG(points(i,2),points(i,1),d+1))) 
            interest_points(count,:) = points(i,:);
            count = count+1;
        end
    end
end

end