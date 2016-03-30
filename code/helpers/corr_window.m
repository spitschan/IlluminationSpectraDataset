function corrSig = correlationWindow(signal, windowSpanUnilateral)
% corrSig = correlationWindow(signal, windowSpanUnilateral);
%
% Calculates a correlation in a given window span
%
% 3/26/2016 	spitschan	Added to repository

if mod(windowSpanUnilateral, 2) ~= 0
   error('Window span must be even.'); 
end

idx1 = [-windowSpanUnilateral:2:-1];
idx2 = [-windowSpanUnilateral+1:2:0];

idx3 = [0:2:windowSpanUnilateral]; idx3(1) = []; % Remove 0
idx4 = [1:2:windowSpanUnilateral];

idxA = [idx1 idx3];
idxB = [idx2 idx4];

%% Calc the correlation
indexVector = 1:length(signal);

for i = 1:length(signal)
   % Find the valid indices
   tmp = all([i+idxA ; i+idxB]' > 1, 2) & all([i+idxA ; i+idxB]' < length(signal), 2);
   corrSig(i) = corr2(signal(i+idxA(tmp)), signal(i+idxB(tmp)));
end

end