function c = calculateC( cmppxx, w, nD )

mmwW = 0.125; % Window width [mm]
AborW = round(mmwW/(cmppxx*10)/2);
c = [];
c(1,1) = 1 + AborW + 4;
c(2,1) = w - AborW - 4;

c1 = c(1);
cE = c(end);
for iD = 1 : nD
    i = iD*2-1;
    c(i) = round( (cE-c1)/(nD-1)*(iD-1)+c1 );
    c(i+1) = c(i);
end