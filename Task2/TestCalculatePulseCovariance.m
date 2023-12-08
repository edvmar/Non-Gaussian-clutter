clear
clc

numberOfPulses = 5;
delta = 1/numberOfPulses;

indeces = 1:numberOfPulses;
jColumn = repmat(indeces, width(indeces),1);
iRow = repmat(indeces',1, width(indeces));

sigma = exp(-(iRow-jColumn).^2*delta)

