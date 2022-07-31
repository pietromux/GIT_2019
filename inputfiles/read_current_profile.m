fid=fopen(output,'r');
tline = fgets(fid);
current=[];
slicenum=[];
ic=1;
is=1;
while ischar(tline)
    tline = fgets(fid);
       if strfind(tline,'current')  
           startIndex=regexpi(tline,'\d');
           current(ic)=str2num(tline(1:startIndex(end)));
               ic = ic+1;
       elseif strfind(tline,'output: slice')
           startIndex=regexpi(tline,'\d');
           slicenum(is)=str2num(tline(startIndex(1):end));
               is = is+1;
       end
end
