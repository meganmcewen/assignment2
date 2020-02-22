%Part 2d
%This code investigates how changing the conductivity of the contacts
%impacts the current density
%Only the conductivity inside the contacts is looked at here - code can
%easily be modified to do this for the area outside the contacts.

clear all

%constants
L = 30; %size of matrix in x
W = L*(2/3); %size of matrix in y
Lb = L/4; %length of contact
Wb = W/4; %width of contact
sigmaOutside = 1; %conductivity outside of contacts; kept constant for now

%create initial matrices
G = sparse(L*W,L*W);
F = zeros(1,L*W);

%change the conductivity inside the contacts each iteration
for sigMult = 1:25 
sigma = ones(W,L);
sigmaInside = sigMult*1e-2;

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
[Ex Ey] = gradient(-VMatrix);
Jx = sigma.*Ex;
Jy = sigma.*Ey;
J_avg = sqrt(Jx.^2 + Jy.^2);
I(sigMult) = sum(J_avg,'all'); %save current for each case of sigma to be plotted

end
figure(1)
plot(1:25,I)
title('Current vs Conductivity Inside Contacts');
xlabel('Sigma Value (x0.01)');
ylabel('Current (A)');
