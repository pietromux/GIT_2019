outputfieldfile = strcat(git_par.dirname,'mg.out.dfl');
Cfieldfile = fopen(outputfieldfile);
Cfield = fread(Cfieldfile,2*git_par.ncar*git_par.ncar*niter,'double');
fclose(Cfieldfile);                                     
iindex = 0;
recfield = Cfield(iindex+1:2:iindex+2*git_par.ncar*git_par.ncar-1);
imcfield = Cfield(iindex+2:2:iindex+2*git_par.ncar*git_par.ncar);
dx = 2*git_par.dgrid/git_par.ncar;
lambda = git_par.lambda;
k = 2*pi/lambda;
propfield = reshape((recfield-1i*imcfield),ncar,ncar); %note that by convention we need the complex coniugate of Genesis field
prop1 = propfield;
clear z_cav
clear rms_out
Ltot=15;
zw=1.4;
Lu=6.304;
f=git_par.f;
d=0.33;
L1=Ltot/2-(Lu-zw)+d;
L2=Ltot/2-zw-d;
atten = 1;
lens = [1 0; -1/f 1]; 
drift1= [1 L1; 0 1];
Mag =1;
ij = 1;
z_cav(1)=Lu;
   power_out= sum(sum(abs(propfield).^2));
   rms_cav(1) = sqrt(sum(sum(abs(propfield).^2.*(x2.^2)))/power_out);
for dz = 0:0.1:L1
   M1 = [1 dz; 0 1];
   field_in_cavity = free_prop_spectr_fun(git_par,M1,Mag,propfield);
   power_out= sum(sum(abs(field_in_cavity).^2));
   ij = ij+1
   z_cav(ij) = Lu+dz;
   rms_cav(ij) = sqrt(sum(sum(abs(field_in_cavity).^2.*(x2.^2)))/power_out);
end
for dz = 0:0.1:L2
   M1 = [ 1 dz; 0 1];
   field_in_cavityb = free_prop_spectr_fun(git_par,M1*lens*drift1,Mag,propfield);
   power_out= sum(sum(abs(field_in_cavityb).^2));
   ij = ij+1
   z_cav(ij) = Lu+L1+dz;
   rms_cav(ij) = sqrt(sum(sum(abs(field_in_cavityb).^2.*(x2.^2)))/power_out);
end
shiftz = Lu+L1-7.5;
z_cav2 = z_cav-shiftz;
z_cav2(z_cav2>7.5) = z_cav2(z_cav2>7.5)-15;
figure(200)
plot(z_cav2,rms_cav*2,'o')
hold on
plot(z+1.6-shiftz,r_size*2,'b')
hold on
plot(1.6-shiftz,init_rms(end)*2*0.9,'bo')

xpl = -7.5:0.1:7.5
w = 0.44e-3*sqrt(1.+(xpl/2.31).^2)
plot(xpl,w,'r')
hold off





