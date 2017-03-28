function S = filterCutFrequency(S,fs,cutFrequency,n)
[z, p, k] = cheby1(n,0.1,2*cutFrequency/fs,'low');
[sos,g]=zp2sos(z,p,k);
h = dfilt.df2sos(sos,g);
S = filter(h,S);
