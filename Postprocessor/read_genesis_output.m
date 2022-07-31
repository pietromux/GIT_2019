lwig = nwig;
nheader = 189+lwig;    % might need to adjust
clear power energy bunching r_size xrms p_midsc z yrms
for js = 1:niter
   C(js) = importdata(output,' ',nheader+(lwig+8)*(js-1));
   power(slippagetot+js,:) = C(js).data(:,1);
   energy(slippagetot+js,:) = C(js).data(:,5);
   bunching(slippagetot+js,:) = C(js).data(:,6);
   p_midsc(slippagetot+js,:) = C(js).data(:,2);
   r_size(slippagetot+js,:) = C(js).data(:,4);
   xrms(slippagetot+js,:) = C(js).data(:,7);
   yrms(slippagetot+js,:) = C(js).data(:,8);
end
