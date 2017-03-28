function [pks,locs,df] = findSpectralPeaks(E,fs,start)

% pick a representative second of noise
e = E(start:(start+fs));
ffe = fft(e);
ffe = ffe(1:fs/2);
df = linspace(0,22050,length(ffe));
in = df<5000;
%x = df(in);
y = log10(abs(ffe(in)));
[pks,locs] = findpeaks(y,'MinPeakDistance',210);