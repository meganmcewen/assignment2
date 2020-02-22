%Part 1b - 2D solution with analytical solution comparison
%This code generates two solutions - one using the G matrix and one using
%the analytical equation. 

clear all

%constants
L = 30; %size of matrix in x
W = 20; %size of matrix in y
iter = 100; %number of iterations
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

        if(x ==1)
            G(n,n) = 1;
            F(n) = 1;
        elseif (x==L)
            G(n,n) = 1;
            F(n) = 1;
        elseif (y ==1)
            G(n,n) = 1;
        elseif(y ==W)
            G(n,n) = 1;
        else
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
title('2D Finite Difference Solution with G matrix');

%analytical solution
a = L;
b = W/2;
x_a = linspace(-b,b,L);
y_a = linspace(0,a,W);

[X,Y] = meshgrid(x_a,y_a);
sum = zeros(size(X));

        for (n = 1:2:iter)
            %use analytical equation to find sum
            sum = sum + ((1/n)*((cosh((n*pi*X)/a))./(cosh((n*pi*b)/a))).*(sin((n*pi*Y)/a)));
            figure(2)
            surf(sum); %watch analytical solution converge
        end
       V_A = ((4*V0)/pi)*sum;
       
figure(3)
surf(V_A); %final analytical solution
title('2D Finite Difference Solution using Analytical Series - Final Convergence');
      




