function newS = noiseRemoveFinal(S,E,fs)

% A multipart noise removal scheme for data formatted according to 
% "Measurement of acoustic and anatomic changes in oral and maxillofacial 
% surgery patients".
% The noise removal is done in several steps: first the speech data is
% adjusted for the frequency response of the speech channel, then  an
% optimized subtraction is done using error channel in  short overlapping 
% time windows, and finally spectral peaks of errors channel are found, and
% these are filtered out using 20th degree Chebyshev Type I filters. 
% 
% Input:    S: the speech data vector
%           E: the MRI-error vector 
%           fs : the audio vector sample rate
%           draw : if you want a fancy graphic
% Output:   newS : the filtered audiodata as single channel scalar signal.
%           fs : the sample rate of newS
            

%% Set parameters
% set the optimization parameters: degree of filter
n = 20; 
% ... and the frequency of final lowpass-filter
cutFrequency = 8000;
% How wide frequency cuts do we want?
cutParameter = 0.005;
%% Linear algebra stuff
E = noisePurifier(S,E,fs);
L = min(length(E),length(S));
S = S(1:L);

%% Filter the frequency response of the waveguide

S = filterFrequencyResponse(S,fs);
%% Start of error signal

start = fs*13;
%for sentences
%start = fs*8;


%% Find spectral peaks
[pks,locs,df] = findSpectralPeaks(E,fs,start);


%% compute the notches, create filter

% If there is only one fundamental frequency use this. It is far less 
% prone to errors. 
%m = floor(mode(fix(diff(df(locs))) ))
% notches = ceil(fs/m);5


m = findPeakHarmonics(pks,locs,df);
%% Filter out the harmonics

S = filterHarmonics(m,S,E,fs,cutParameter);


%% Take out freqs above wanted

S = filterCutFrequency(S,fs,cutFrequency,n);
% if necessary - usually not, and kind of questionable, but efficient for
% sure
S = finalizeFiltering(S,E,fs);


%% Finally (if necessary) take out stationary noise
newS = spectralSubtraction(S,fs,6);

