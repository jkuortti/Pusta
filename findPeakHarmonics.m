function m = findPeakHarmonics(pmags,locs,f)
% m will be the vector containing the notch numbers/harmonic
widthParameter=7;


m=[];
plocs = f(locs);


while ~isempty(plocs)

[v,c] = max(pmags);
p = plocs(c);
plocs(c) = [];
pmags(c) = [];
[Q,V] = sort(abs(p-plocs),'ascend');

ptemp = plocs(V);
D = abs(ptemp - p);
ptemp(D<90)= [];
D = abs(ptemp - p);
foundHarmonics = [];
for k = 1:length(ptemp)
    foundTemp =  find(min(mod(fix(D'),fix(D(k))), abs(mod(fix(D'),fix(D(k)))-D(k)) ) < widthParameter);
    if length(foundTemp)>3
        foundHarmonics = union(foundHarmonics,foundTemp);
    end
end

% There WILL be variance in foundHarmonics; hoping median will be robust
% enough to identify "the center".
m = [m,median(diff(sort(ptemp(foundHarmonics))))];%#ok - will not grow much
plocs = setdiff(plocs,ptemp(foundHarmonics));
[plocs,indL] = setdiff(plocs,ptemp(foundHarmonics));
pmags = pmags(indL);
end