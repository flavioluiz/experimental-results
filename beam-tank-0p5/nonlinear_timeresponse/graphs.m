clear all;
load symexc_5thmode_3V
figure;
plot(outputvalues(:,1)/3); hold all;
load symexc_5thmode_5V
plot(outputvalues(:,1)/5);

load symexc_5thmode_10V
plot(outputvalues(:,1)/10);

load symexc_5thmode_15V
plot(outputvalues(:,1)/15);

legend('3 V', '5 V', '10 V', '15 V');