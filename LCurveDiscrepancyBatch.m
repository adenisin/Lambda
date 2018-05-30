function [lambda,lambda_Discrepancy,files] = LCurveDiscrepancyBatch(E,conversion,NA, PSFwidth, Fmin)
%calculates tfm between 2 frames
nu = 0.49;
A=(1+nu)/(pi*E);
searchitem = input('Type the base of your files (i.e., *20kPa*): ');
files = dir(searchitem);
%E = 18300; 46700 %Pa stiffness
%nu = 0.5; %Poission's ratio
%conversion = 0.163; %um/pixel

lambda = zeros(size(files,1),1);
lambda_Discrepancy = zeros(size(files,1),1);
lambda_sq = zeros(size(files,1),1);

for z = 1:size(files,1)
    data = dlmread(files(z).name);
    x = data(:,1);
    y = data(:,2);
    u = data(:,3)*conversion*1e-6;
    v = data(:,4)*conversion*1e-6;
    mag = data(:,5)*conversion*1e-6;
    
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
    
    %%determine lambda
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
    [U,s,~]=csvd(Gt);
    
    [lambda_sq(z),rho(:,:,z),eta(:,:,z),reg_param_sq(:,z)]=l_curve(U,s,ut(end:-1:1));
    [reg_c,rho_c,eta_c] = l_corner(rho(:,:,z),eta(:,:,z),reg_param_sq(:,z),U,s,ut(end:-1:1));
    lambda(z) = sqrt(lambda_sq(z));
    reg_param(:,:,z) = sqrt(reg_param_sq(:,z));
    
    % %need to apply the Morozov Discrepancy Priciple to find a lambda which may be
    % %better than lambda at the corner for our data
    % %based on the chi sq distribution with specific degrees of freedom
    % %following Schwarz 2002
    %need to apply the Morozov Discrepancy Priciple to find a lambda which may be
    %better than lambda at the corner for our data
    %based on the chi sq distribution with specific degrees of freedom
    %following Schwarz 2002
    
    %sigma = PSF_spread*conversion;
    %Schwarz 2002 defines sigma (standard deviation) as sqrt (4(N-M))
    %Assigning degrees of freedom
    %N is the number of sites at which there are focal adhesions present
    %M is the positions at which force can be measured
    N = max(length(v), length(u)); %displacement field measured at different sites
    %M is number of sites which are actively making contractions
    %chi squard = mag(GF = u)^2/sigma^2, where sigma is standard deviation of
    %measurement errors for displacement u
    %choose residual norm such that R = 2(N-M)*sigma^2
    %set minimum force
   % Fmin = 2.5E-9; %nN, 1E-9 kgm/s^2, focal adhesions are 100 nN
    %look at magnitude vector (mag) to detremine how many of the
    %displacement vectors are signal versus noise
    %use PSF for our microscope to get minimum resolution
    R = (548E-9)/(2*NA); %Abbe's diffraction limit %RFP: 548E-9 m, GFP: 488E-9 m
    d = Fmin*(1+nu)/(R*E*pi);
    signal = mag>d;
    sumsignal = sum(signal);
    M = (sumsignal/length(mag))*N;
    %sigma is the standard deviation of measurement errors for displacement
    %define sigma as the width of the PSF, gaussian fit to a single particle
    sigma = PSFwidth*conversion;
    rho_opt_abs = 2*(N-M)*(sigma^2);
    rho_opt = 2*(N-M)*(sigma^2)*(conversion*1E-6)^2;
    
    %recreate the L curve using Tikhonov
    lambda_range =[linspace(1E-15, 1E-13, 200),linspace(1E-13, 1E-11, 200),linspace(1E-11, 1E-8, 200)]';
    [lambda_f,rho_range,eta_range] = tikhonov (U,s,V,ut(end:-1:1),(lambda_range).^2);
    
    val = abs(rho_range-rho_opt);
    [idx idx] = min(val); %index of closest value
    lambda_Discrepancy(z) = lambda_range(idx);
end
disp('finished');
end

