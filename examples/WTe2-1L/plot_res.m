% hr = HR.from_wannier90('wannier90_symmed_hr.dat');
% hr.bandplot();
% 
% hr_Bxy = add_zeeman_field(hr, [200 200 0]);
% hr_Bxy.Gen_hr('wannier90_symmed_hr_Bxy.dat');
% hr_Bxy.bandplot();
%%
df = readmatrix('sigma_SOAHC_int_eta1.00meV.dat','NumHeaderLines',3);
figure()
hold on
plot(df(:,1), df(:,2),'LineWidth',1.5);
plot(df(:,1), df(:,3),'LineWidth',1.5);
hold off

xlabel('E(eV)')
xlim([-0.03 0.05])
ylabel("\chi")
legend(["\chi_{xyy}", "\chi_{yxx}"])
set(gca,'FontSize',21)
