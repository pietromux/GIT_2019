read_current_profile;
read_genesis_output;
p_out = power;
figure(30) 
inum = 1:size(slicenum,2);
plot(inum+slippagetot,current)
hold on
plot(p_out(:,end)/1e6)
hold off
%%
z = (0:nwig)*lambda_w;
K_undulator = [lat(1,:) 0];

figure(1)
plot(z,power)
title('power')

figure(2)
plot(z,energy)
title('Average beam energy')

figure(4)
plot(z,K_undulator(1:length(z)))
title('K undulator')

figure(6)
plot(z,bunching)
title('bunching')

figure(7)
plot(z,p_midsc)
title('On axis intensity')

figure(3)
plot(z,r_size,'g')
hold on
plot(z,xrms,'r')
plot(z,yrms,'k')
hold off
title('Beam and radiation rms spot size')

% outres.z= z;
% outres.K_undulator = K_undulator;
% outres.lambda_undulator = lambda_undulator;
% outres.power = power;
% outres.phislice = phislice;
% outres.energy = energyslice;
% outres.bunching = bunching;
% outres.trapped = ncslice;
% outres.p_mid = p_mid;
% outres.r_size = r_size;
% 
% outres.xrms = xrms;
% outres.d3fac = d3fac;
% outres.resphase = resphase;
% outres.gammaresonant = gammaresonant;
% outres.theta = thetaslice;
% outres.deltatheta = deltatheta;


