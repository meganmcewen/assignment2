%Part 1a - 1D solution 
%This code is the simplest case for the finite difference code, looking at a 1D solution 

clear all

%constants
L = 30; %size of matrix in x
W = 20; %size of matrix in y
V0 = 1; %initial velocity

%create initial matrices
G = sparse(L*W,L*W);
F = zeros(L*W,1);

%generate G matrix
for x = 1:L
    for y = 1:W
        
    %mapping equation
    n = y + (x-1) * W;

    %local mapping
    nxm = y+(x-2)*W;
    nxp = y+(x)*W;
    nym = (y-1)+(x-1)*W;
    nyp = (y+1)+(x-1)*W;

        if(x ==1) %left x boundary
            G(n,n) = 1;
            F(n) = V0;
        elseif (x==L) %right x boundary
            G(n,n) = 1;
            F(n) = 0;
        elseif (y ==1) %bottom y boundary
           G(n,:) = 0;
           G(n,nxm) = 1;
           G(n,n) = -3;
           G(n,nxp) = 1;
           G(n,nyp) = 1;
        elseif(y ==W) %top y boundary
           G(n,:) = 0;
           G(n,nxm) = 1;
           G(n,n) = -3;
           G(n,nxp) = 1;
           G(n,nyp) = 1;
        else %all interior points
           G(n,:)= 0;
           G(n,nxm) = 1;
           G(n,n) = -4;
           G(n,nxp) = 1;
           G(n,nym) = 1;
           G(n,nyp) = 1;
        end
    end
end

V = G\F; %use G matrix to solve for V

%populate V matrix solution
for (x = 1:L)
    for (y = 1:W)
        n = y + (x-1) * W;
        VMatrix(x,y) = V(n);
    end
end

%plot g matrix solution
figure(1)
surf(VMatrix);
title('1D Finite Difference Solution with G matrix');

