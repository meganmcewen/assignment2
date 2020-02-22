%Part 2b
%This code observes what happens when the mesh size is varied, but
%conductivity and bottleneck size remain consistent.

clear all

for meshMulti = 1:5 
%can do this for more meshes - kept this to a small number for testing 
%due to the time it takes for the solution to run.
    
%constants
L = 30*meshMulti; %size of matrix in x
W = L*(2/3)*meshMulti; %size of matrix in y
Lb = L/4; %length of contact
Wb = W/4; %width of contact
sigmaInside = 1e-2; %sigma value for contacts
sigmaOutside = 1;

%create initial matrices
G = sparse(L*W,L*W);
F = zeros(1,L*W);
sigma = ones(W,L);
    
%create contacts
contactLeft = (L/2) - (Lb/2);
contactRight = (L/2) + (Lb/2);
contactBottom = Wb;
contactTop = W - Wb;

%set up sigma for contacts
sigma(:,:) = sigmaOutside;
sigma(contactTop:W, contactLeft:contactRight) = sigmaInside;
sigma(1:contactBottom, contactLeft:contactRight) = sigmaInside;

%G matrix
for x=1:L
    for y=1:W
        
        %mapping equation
        n = y +(x-1)*W;
        
        if(x==1) 
         
            G(n,n) = 1;
            F(n)= 1;
            
        elseif(x==L)
            
            G(n,n) = 1;
            F(n)=0;
            
        elseif(y == 1)
            
            %local mapping
            nyp = y+1+(x-1)*W;
            nxp = y+(x)*W;
            nxm = y+(x-2)*W;

            sig_yp =(sigma(y,x)+sigma(y+1,x))/2;
            sig_xp=(sigma(y,x)+ sigma(y,x+1))/2;
            sig_xm =(sigma(y,x)+ sigma(y,x-1))/2;
            
            G(n,n)= -(sig_yp+sig_xp+sig_xm);
            G(n,nyp)= sig_yp;
            G(n,nxp)=sig_xp;
            G(n,nxm)= sig_xm;
            
        elseif(y==W)
            
            %local mapping
            nxp = y+(x)*W;
            nxm = y+(x-2)*W;
            nym = y-1+(x-1)*W;

            sig_xp=(sigma(y,x)+ sigma(y,x+1))/2;
            sig_xm =(sigma(y,x)+ sigma(y,x-1))/2;
            sig_ym =(sigma(y,x)+ sigma(y-1,x))/2;
           
            G(n,n)=-(sig_ym+sig_xp+sig_xm);
            G(n,nym)=sig_ym;
            G(n,nxp)=sig_xp;
            G(n,nxm)=sig_xm;
            
        else
            
            %local mapping
            nyp = y+1+(x-1)*W;
            nxp = y+(x)*W;
            nxm = y+(x-2)*W;
            nym = y-1+(x-1)*W;

            sig_yp =(sigma(y,x)+sigma(y+1,x))/2;
            sig_xp=(sigma(y,x)+ sigma(y,x+1))/2;
            sig_xm =(sigma(y,x)+ sigma(y,x-1))/2;
            sig_ym =(sigma(y,x)+ sigma(y-1,x))/2;
        
            G(n,n)=-(sig_yp+sig_ym+sig_xp+sig_xm);
            G(n,nyp)= sig_yp;
            G(n,nym)= sig_ym;
            G(n,nxp)= sig_xp;
            G(n,nxm)= sig_xm;
            
        end
    end
end

V = G\F';

VMatrix = zeros(W,L);

%populate V matrix solution
for (x = 1:L)
    for (y = 1:W)
        n = y + (x-1)*W;
        VMatrix(y,x) = V(n);
    end
end

%solve for J
[Ex,Ey] = gradient(-VMatrix);
Jx = sigma.*Ex;
Jy = sigma.*Ey;
J_avg = sqrt(Jx.^2 + Jy.^2);
I(meshMulti) = sum(J_avg,'all'); %current for each mesh is saved to be plotted

end

%plots
figure(1)
plot(1:5,I) %change this based on number of meshes testing for
hold on
title('Mesh Density and Current')
xlabel('Mesh Multiplier')
ylabel('Average Current (A)')
grid on