%Part 2c
%This code looks at how narrowing the bottleneck impacts the current
%density
 
clear all

%constants
L = 30; %size of matrix in x
W = L*(2/3); %size of matrix in y
Lb = L/4; %length of contact
Wb = W/4; %width of contact
sigmaInside = 1e-2; %sigma value for contacts
sigmaOutside = 1;

%create initial matrices
G = sparse(L*W,L*W);
F = zeros(1,L*W);
sigma = ones(W,L);

%bottleneck varying
for narrow = 1:10 %testing 10 iterations of increasingly narrow bottleneck
    
%create contacts
contactLeft = (L/2) - (Lb/2);
contactRight = (L/2) + (Lb/2);
contactBottom = Wb + (narrow*.5);
contactTop = (W - Wb)-(narrow*.5);

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
I(narrow) = sum(J_avg,'all'); %save current for each bottleneck to plot at end

figure(1)
surf(sigma); %plot this to watch the bottleneck narrow
view(2)

bottleNeck(narrow) = (contactTop - contactBottom); %save bottleneck sizes

end

%plots
figure(2)
plot(bottleNeck,I)
title('Narrowing Bottleneck and Current')
xlabel('Bottleneck Size')
ylabel('Average Current (A)')
grid on