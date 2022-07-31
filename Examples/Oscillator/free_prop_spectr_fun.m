function [ propfield ] = free_prop_spectr( git_par, M, Mag, field )

%% Spectral fourier Wavefront Propagation algorithm
% Inputs 
% field  = Complex input field in space domain
% dx space step of the field array
% M propagation matrix
% Mag magnification of final grid
%%%%%%%
% Output
% propfield = complex output field in space domain
% only works for single slice

dx = 2*git_par.dgrid/git_par.ncar;
lambda = git_par.lambda;
k = 2*pi/lambda;
 
    A = M(1,1);
    B = M(1,2);
    C = M(2,1);
    D = M(2,2);

    [Mx, My] = size(field);
    a1 = dx*(Mx+1)/2;
    a3 = a1;
    a2 = Mag*a1;
    a4 = a2;

    dkx = 2*pi/(k*dx*Mx) ; dky = 2*pi/(k*dx*My); 
    nx = ((1:Mx)-(Mx+1)/2);
    ny = ((1:My)-(My+1)/2);
    [kx,ky] = meshgrid(nx*dkx,ny*dky);
    [x1,y1] = meshgrid(nx*dx,ny*dx);
    [x2,y2] = meshgrid(nx*dx*Mag,ny*dx*Mag);

    vfield = field.*exp(-1i*(pi*(A-Mag).*(x1.^2+y1.^2)/B/lambda));
    Fresnel = Mag/B/lambda;
    psift = (fftshift(fft2(vfield))); 
    psi_k = psift.*exp(1i.*(kx.^2+ky.^2)*k^2/4/pi/Fresnel);
    propvfield = (ifft2(ifftshift(psi_k)));
    propfield = propvfield.*exp(-1i*(pi*(D-1/Mag).*(x2.^2+y2.^2)/B/lambda));
    %% output
%    power_in = sum(sum(abs(field).^2));
%    power_out= sum(sum(abs(propfield).^2));
%    rms_in = sqrt(sum(sum(abs(field).^2.*(x1.^2+y1.^2)))/power_in);
%    rms_out = sqrt(sum(sum(abs(propfield).^2.*(x2.^2+y2.^2)))/power_out);
%    formatSpec = 'Rad_propagation. Power %.3e ->  %.3e.   \n rms %.3e -> %.3e . \n';  
%    fprintf(formatSpec,power_in,power_out,rms_in,rms_out);
end
    
    
