function filteredS = filterFrequencyResponse(S,fs)
% load magnitudal frequency responses from file, construct fir-filter and
% filter the signal S, then return.


load('puhevasteMag.mat');
puhevaste = decimate(puhe_response,80);
f = downsample(puhe_f,80);
puhevaste = [puhevaste(1) puhevaste puhevaste(end) puhevaste(end)];
f = [0  f 4500 22050];
%[b, a] = fir2(90,2*f./fs,-puhevaste); %doesn't work, amirite
[b,a] = yulewalk(8,2*f./fs,-puhevaste); 
filteredS = filter(b,a,S);