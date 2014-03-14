a=[0.35 0.8 1.2 1.2 0.8 0.35 0.35]
b=[2.5 2.5 5.5 7.5 7.5 5 2.5]

plot(a,b,0.867,4.45,"o");
xlabel('omega nsp');
ylabel('zeta sp (rad/s)');
xlim([0,1.4]);
ylim([0,8]);

grid on;

print("output/test.png",'-dpng','-r100')