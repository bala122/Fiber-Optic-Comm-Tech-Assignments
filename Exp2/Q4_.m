%Q4
%Generating 2^12 prbs symbols using the functions prbs of order 9.
%Directly done because the prbs function repeats the symbols after 2^O -1
%symbols.
no_symb_orig = 2^12;
y_prbs9= prbs(9,no_symb_orig);
scatterplot(y_prbs9)
grid on
xlabel('I-phase')
ylabel('Q-phase')
title('Constellation diagram for OOK(no noise)')
legend('OOK (no noise)')



%Q4
%Using the Nt = 8 signal for a better prbs waveform
%converting  input signals to matrices
sig5=zeros(length(c3),1);
sig5 = transpose(c3);
constdiag = comm.ConstellationDiagram('Title','Constellation diagram for OOK(no noise)');
constdiag(sig5)

%Q5 adding noise
%5dB
out_5dB = awgn(sig5, 5);
%scatterplot(out_5dB)
constdiag2 = comm.ConstellationDiagram('Title','Constellation diagram for OOK(SNR 5dB)');
constdiag2(out_5dB)

%10dB
out_10dB = awgn(sig5,10);
constdiag3 = comm.ConstellationDiagram('Title','Constellation diagram for OOK(SNR 10dB)');
constdiag3(out_10dB)


%15dB

out_15dB = awgn(sig5,15);
constdiag4 = comm.ConstellationDiagram('Title','Constellation diagram for OOK(SNR 15dB)');
constdiag4(out_15dB)


