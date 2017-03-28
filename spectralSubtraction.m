function y = spectralSubtraction(x,fs,IS)

% Implements the spectral subtraction algorithm from the article of S.
% Boll  Suppression of acoustic noise in speech using spectral subtraction 

% Input arguments:
% x : the original signal
% fs : the sampling rate of the original signal
% IS : the length of the initial silence in seconds (noise only)

% Inspired by the codes of 
% Esfandiar Zavarehei
% April -05
% https://se.mathworks.com/matlabcentral/fileexchange/7675-boll-spectral-subtraction
% and
% Philipos Loizou
% Speech Enhancement: Theory and Practice

wlen = fix(.025*fs); %Window length is 25 ms
% Shift percent is 40% (10ms)
overlapPercent = .5; 

% number of initial silence segments
NIS = fix((IS*fs-wlen)/(overlapPercent*wlen) +1);

% divide signal into: timeframes-> output matrix y has frames in
% columns and compute spectrogram of the signal
Y = fft(segment(x,wlen,overlapPercent),wlen);
theta = angle(Y(1:fix(end/2)+1,:)); 
Y = abs(Y(1:fix(end/2)+1,:)); 

NFrames = size(Y,2);
NoiseEstimate = mean(Y(:,1:NIS),2); %initial Noise Power Spectrum mean
NRM = zeros(size(NoiseEstimate));% Noise Residual Maximum (Initialization)
NoiseCounter = 0;

MagnitudeAverage = movmean(Y,3,2);


X = zeros(size(Y));
for k=1:NFrames
    [isSpeech, NoiseCounter]=vad(Y(:,k),NoiseEstimate,NoiseCounter); 
    if isSpeech
        D = MagnitudeAverage(:,k) - NoiseEstimate; % Specral Subtraction
        if k>1 && k<NFrames
            aux = min([D MagnitudeAverage(:,k-1)-NoiseEstimate MagnitudeAverage(:,k+1)-NoiseEstimate],[],2);
            D(D<NRM) = aux(D<NRM);
        end
        X(:,k) = max(D,0);
    else
        % this window is noise
        NoiseEstimate = 0.9*NoiseEstimate+0.1*Y(:,k);% Update and smooth noise
        NRM = max(NRM,MagnitudeAverage(:,k)-NoiseEstimate);% Update Maximum Noise Residue
        SNRseg = 20*log10(norm(Y(:,k))/norm(NoiseEstimate)); 
        Beta = berouti1(SNRseg);
        X(:,k) = Beta*Y(:,k);
    end
end
y = OverlapAdd(X,theta,wlen,overlapPercent*wlen);
y = y./max(y);



function Y = OverlapAdd(X,theta,winLen,shiftLen)
% Y is the time domain signal reconstructed from spectrogram X
% X is a (magnitude) spectrogram of a signal with shiftLen overlap
% theta is the phase angle of the spectrum which should have the same dimension as X. 
% winLen is the window length of time domain segments in samples
% shiftLen is the overlap fraction of the frames in X. 

NFrames=size(X,2);
Spectrum = X.*exp(1i*theta);
isEven = rem(winLen,2);
Spectrum = [Spectrum;flipud(conj(Spectrum(2:end-1+isEven,:)))];
shiftLen=fix(shiftLen);
Y = zeros((NFrames-1)*shiftLen+winLen,1);
for k=1:NFrames
    start = (k-1)*shiftLen+1;
    Y(start:start+winLen-1) = Y(start:start+winLen-1)+real(ifft(Spectrum(:,k),winLen));
end



function [isSpeech,NoiseCounter] = vad(signal,noise,NoiseCounter)
% Spectral Distance Voice Activity Detector
% Signal : is the the current frames magnitude spectrum which is to labeled as
%          noise or speech, 
% Noise :  is noise magnitude spectrum template (estimation),
% NoiseCounter : counts the number of consecutive frames of pure noise.  
% isSpeech : is True if frame contains speech, False if not.  

NoiseMargin = 3; NoisePeriodAssumption = 6;
sdist= mean(max(20*(log10(signal)-log10(noise)),0));
NoiseCounter = (sdist<NoiseMargin)*(NoiseCounter+1);
isSpeech = (NoiseCounter<NoisePeriodAssumption);



function S = segment(signal,winlen,overlap)
% SEGMENT chops a signal to overlapping windowed segments
% A= SEGMENT(X,W,SP,WIN) returns a matrix which its columns are segmented
% and windowed frames of the input one dimentional signal, X.
% W is the % number of samples per window. 
% overlap is the overlap percentage,

Window = hamming(winlen); 
overlap = fix(winlen.*overlap);
N = fix((length(signal)-winlen)/overlap +1); %number of segments
Index = (1:winlen)' + ((0:(N-1))*overlap);
S = bsxfun(@times,signal(Index),Window);



function a = berouti1(SNR)
% From Speech Enhancement 
% by Philipos Loizou
if SNR>=-5.0 && SNR<=20
   a=3-SNR*2/20;
else
  if SNR<-5.0
   a=4;
  end
  if SNR>20
    a=1;
  end
end