function S = finalizeFiltering(S,E,fs)

ffs = fft(E(fs*13:fs*14));
ffs = ffs(1:fs/2);
x = linspace(0,fs/2,length(ffs));
y = 10*log10(abs(ffs));
[~,locs] = findpeaks(y,'MinPeakDistance',210);
frqs = x(locs);

for k = 1:length(frqs)
    w0 = frqs(k)/(fs/2);  bw = w0/60;
    [b,a] = iirnotch(w0,bw);
    S = filter(b,a,S);
end