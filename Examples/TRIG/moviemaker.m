%MovieTRIG(1)=[];
myVideo = VideoWriter('TRIG_ps.avi');
myVideo.FrameRate = 5;
open(myVideo);
writeVideo(myVideo, MovieTRIG);
close(myVideo);