clear;
clc;

syms P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13 P14 P15

P = [P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15];

DeltaN = [[-1, 1, zeros(1, 13)]', [0, -1, 1, zeros(1, 12)]', \
     [zeros(1, 2), -1, 1, zeros(1, 11)]', [zeros(1, 3), -1, 1, 1, zeros(1, 9)]',\ 
     [zeros(1, 4), -1, 0, 1, zeros(1, 8)]', [zeros(1, 5), -1, 1, zeros(1, 8)]',\ 
     [zeros(1, 6), -1, 1, zeros(1, 7)]', [zeros(1, 7), -1, 1, zeros(1, 6)]',\ 
     [zeros(1, 8), -1, 1, 1, zeros(1, 4)]', [zeros(1, 9), -1, 0, 1, zeros(1, 3)]',\
     [zeros(1, 10), -1, 1, zeros(1, 3)]', [zeros(1, 11), -1, 1, 1, 0]',\ 
     [zeros(1, 12), -1, 0, 1]', [zeros(1, 13), -1, 1]', zeros(1, 15)'];

mut = DeltaN*P';

sigmat = zeros(15, 15);
for i = 1:15
  sigmat = sigmat + P(i)*DeltaN(1:15, i)*DeltaN(1:15, i)';
endfor