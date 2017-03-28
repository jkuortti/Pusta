function newE = noisePurifier(S,E,fs)
% Attempt clean speech-and speech like components from 
% error signal E. 
% Input: 	S - signal containing speech, or at least mostly 
%			speech
%			E - Signal to be cleaned, ie. the error signal
%			fs - sample rate for both E and S.
% Output: 	newE - cleaned error signal
w = 0.1*fs;
s = newSegment(S,w);
e = newSegment(E,w);
N = size(e,2);
for k = 1:N
   ek = e(:,k);
   sk = s(:,k);
   ek = ek - (dot(sk,ek/norm(ek)))*ek;
   e(:,k) = ek;
end
newE = e(:);
newE = newE(1:length(S));

