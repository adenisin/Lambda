function [lambda,lambda_Discrepancy] = LCurveDiscrepancy(E,conversion)
nu = 0.49;
A=(1+nu)/(pi*E);
filename  = input('Type the name of your file: ');
%E = 18300; 46700 %Pa stiffness
%nu = 0.5; %Poission's ratio
%conversion = 0.163; %um/pixel
data = dlmread(filename);
x = data(:,1);
y = data(:,2);
u = data(:,3)*conversion*1e-6;
v = data(:,4)*conversion*1e-6;
mag = data(:,5);

%For u and v, need to make a matrix
spacing = (x(2,1)-x(1,1));
u2 = zeros((max(x)-x(1,1))/spacing+1,(max(y)-y(1,1))/spacing+1);
v2 = zeros((max(x)-x(1,1))/spacing+1,(max(y)-y(1,1))/spacing+1);

k = 1;
%this is to take our datapoints from columns into spatial distributions
for j = 1:size(u2,2)
    for i = 1:size(u2,1)
        u2(i,j) = u(k);
        v2(i,j) = v(k);
        k=k+1;
    end
end

%restrict to 50x50 area for performance
nx=min(50,size(u2,1));
ny=min(50,size(u2,2));

%fourier transform
u=(fft2(u2(1:nx,1:ny)));
v=(fft2(v2(1:nx,1:ny)));

%initialize traction vectors w. nans
Tx=zeros(size(u));
Ty=zeros(size(v));

%image size
Nx=size(u,1);
Ny=size(v,2);

%distance between grid points
D=(x(2,1)-x(1,1))*conversion*1e-6; %in m

%wavenumbers: theory from fourier transform; 
%matlab starts at f=0:N/2 then from -N/2 to -1;
dkx=1/(Nx*D);
kx=[0:fix(Nx/2),-fix(Nx/2):-1]*dkx*2*pi;
dky=1/(Ny*D);
ky=[0:fix(Ny/2),-fix(Ny/2):-1]*dky*2*pi;

%determine lambda for regularization
ii=0;
n=0;
Gt=zeros(2*(length(kx)-1)*(length(ky)-1),2*(length(kx)-1)*(length(ky)-1));
ut=zeros(2*(length(kx)-1)*(length(ky)-1),1);
%loop over wavenumbers
for i=kx(1:end-1)
    ii=ii+1;
    jj=0;
    for j=ky(1:end-1)
        jj=jj+1;
        n=n+1;
        
        k=sqrt(i^2+j^2);
        un=u(ii,jj);
        vn=v(ii,jj);
        uk=[un;vn];
        
        %at f=0 the tractions must be 0
        if (i==0) && (j==0)
            Tx(ii,jj)=0;
            Ty(ii,jj)=0;
            
        else
            %Qingzong: at the nyquist frequency, set the off-diagonal element to zero, see Butler et al. Am J Physil Cell Physiol 2001
            if (ii==Nx/2)||(jj==Ny/2)
                K=A*2*pi/(k^3)*[(1-nu)*k^2+nu*j^2,0; 0,(1-nu)*k^2+nu*i^2];
            else
                K=A*2*pi/(k^3)*[(1-nu)*k^2+nu*j^2,-nu*i*j; -nu*i*j,(1-nu)*k^2+nu*i^2];
            end
            
            %regul. for l-curve
            Gt(2*n-1:2*n,2*n-1:2*n)=K'*K;
            ut(2*n-1:2*n,1)=K'*uk;
        end
    end
end

%decompose into singular value matrix and use regutools to solve for
%optimum lambda
disp('Decomposing into SVD...');
[U,s,V] = csvd(Gt);
disp('Decomposed into SVD');
disp(' ');

% use L-curve
[lambda_sq,rho,eta,reg_param_sq]=l_curve(U,s,ut(end:-1:1));
[reg_c,rho_c,eta_c] = l_corner(rho,eta,reg_param_sq,U,s,ut(end:-1:1));
lambda = sqrt(lambda_sq);
reg_param = sqrt(reg_param_sq);

figure(1)
figure1 = loglog(rho,eta);
xlabel('residual norm || A x - b ||_2');
ylabel('solution norm || x ||_2');
figureHandle = gcf;
% make all text in the figure to size 15
set(findall(figureHandle,'type','text'),'fontSize',15)
hold on;
figure1 = plot(rho_c,eta_c, '*r');

figure(2)
figure2 = loglog(reg_param, rho);
xlabel('Lambda');
ylabel('residual norm || A x - b ||_2');
hold on;
figure2 = plot(lambda,rho_c, '*r');
figureHandle = gcf;
%# make all text in the figure to size 20
set(findall(figureHandle,'type','text'),'fontSize',15)
% 
% 
%need to apply the Morozov Discrepancy Priciple to find a lambda which may be
%better than lambda at the corner for our data
%based on the chi sq distribution with specific degrees of freedom
%following Schwarz 2002
%rankG = rank(Gt);
DF = size(Gt,1)-1;
%PSF_spread = 2.25;
%sigma = PSF_spread*conversion;
sigma = 2*conversion;
rho_opt = sqrt(DF)*sigma;
rho_opt_plus = (sqrt(DF)+(2*DF)^0.25)*sigma;
rho_opt_minus = (sqrt(DF)-(2*DF)^0.25)*sigma;

%absolute value of rho
rho_abs = rho/((conversion*1e-6)^2);
rho_c_abs = rho_c/((conversion*1e-6)^2);
figure(3)
figure3 = loglog(reg_param, rho_abs);
xlabel('Lambda');
ylabel('residual norm || A x - b ||_2');
hold on;
figure3 = plot(lambda,rho_c_abs, '*g');
figureHandle = gcf;
%# make all text in the figure to size 20
set(findall(figureHandle,'type','text'),'fontSize',15)
figure3 = plot(reg_param,rho_opt*ones(size(reg_param)), 'r');
figure3 = plot(reg_param,rho_opt_plus*ones(size(reg_param)), 'r:');
figure3 = plot(reg_param,rho_opt_minus*ones(size(reg_param)), 'r:');

%find the lambda at the rho_opt
    val = abs(rho_abs-rho_opt); 
    [idx idx] = min(val); %index of closest value
    lambda_Discrepancy = reg_param(idx);
    
% % preview force measurements
% %fourier transform
% u=(fft2(u2));
% v=(fft2(v2));
% 
% %initialize traction vectors w. nans
% Tx=zeros(size(u));
% Ty=zeros(size(v));
% 
% %image size
% Nx=size(u,1);
% Ny=size(v,2);
% 
% %wavenumbers: theory from fourier transform; matlab starts at f=0:N/2 then
% %from -N/2 to -1;
% dkx=1/(Nx*D);
% kx=[0:fix(Nx/2),-fix(Nx/2):-1]*dkx*2*pi;
% dky=1/(Ny*D);
% ky=[0:fix(Ny/2),-fix(Ny/2):-1]*dky*2*pi;
% 
% ii=0;
% %loop over wavenumbers
% for i=kx(1:end-1)
%     ii=ii+1;
%     jj=0;
%     for j=ky(1:end-1)
%         jj=jj+1;
%         
%         k=sqrt(i^2+j^2);
%         un=u(ii,jj);
%         vn=v(ii,jj);
%         uk=[un;vn];
%         
%         %at f=0 the tractions must be 0
%         if (i==0) && (j==0)
%             Tx(ii,jj)=0;
%             Ty(ii,jj)=0;
%             
%         else
%             %Qingzong: at the nyquist frequency, set the off-diagonal element to zero, see Butler et al. Am J Physil Cell Physiol 2001
%             if (ii==Nx/2)||(jj==Ny/2)
%                 K=A*2*pi/(k^3)*[(1-nu)*k^2+nu*j^2,0; 0,(1-nu)*k^2+nu*i^2];
%             else
%                 K=A*2*pi/(k^3)*[(1-nu)*k^2+nu*j^2,-nu*i*j; -nu*i*j,(1-nu)*k^2+nu*i^2];
%             end
%                        
%             %2D identity
%             H=eye(2);
%             
%             %now finally, the traction force calculation
%             Gn=K'*K+lambda^2*H;
%             T=Gn\(K'*uk(end:-1:1));
%             Tx(ii,jj)=T(2);
%             Ty(ii,jj)=T(1);
%         end
%     end
% end
% 
% %transform back into real space
% Trx=real((ifft2(Tx)));
% Try=real((ifft2(Ty)));
% 
% %norm of the tractions
% v = sqrt( Trx.^2 + Try.^2 );
% 
% %Interpolate traction data for every pixel
% try
% [Xqu,Yqu]=meshgrid(x1(1):x1(end),y1(1):y1(end));
% V = interp2(x1,y1,v,Xqu,Yqu,'linear');
% catch
%     V=v;
% end
% figure, imagesc(V);

end

