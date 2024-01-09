clear
clc

n = 50000;
ROC100 = 'Theoretical best ROC curve';
ROC50  = 'ROC curve for a random 50/50 detector';
ROC0   = 'Theoretical worst ROC curve';

xValues100 = [0.01];
yValues100 = [0.99];

xValues50  = linspace(0,1,n);
yValues50  = xValues50;

xValues60  = linspace(0,1,n);
yValues60  = log10(8.8*xValues60 + 1);

xValues70  = linspace(0,1,n);
yValues70  = log(980*xValues60 + 1)/log(1000);

xValues0   = [linspace(0,1,n), 1];
yValues0   = [zeros(1,n) 1];



figure(1)
hold on
p1 = plot(xValues100, yValues100, 'o');
set(p1, 'markerfacecolor', get(p1, 'color'));

p2 = plot(xValues50, yValues50,'--');
p3 = plot(xValues60, yValues60,'-');
p4 = plot(xValues70, yValues70,'-');
axis([-0.005 1.005 -0.005 1.005])
xlabel('P_{FA}')
ylabel('P_{TD}')
%legend(ROC100, ROC50)
hold off



