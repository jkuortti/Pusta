function ns = filterHarmonics(m,S,E,fs,cutParameter)
ns = S;
es = E;

for k = 1:length(m)
    notches = ceil(fs/m);
    D = fdesign.comb('notch','N,BW',notches,cutParameter);
    H  = design(D);
    a = H.Denominator;
    b = H.Numerator;
    ns = filter(b,a,ns);
    es = filter(b,a,es);
end
