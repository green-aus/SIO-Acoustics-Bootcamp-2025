%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Bootcamp Record
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%% Wave file names
fname1='test2_03.wav';
fname2='test2_002.wav';
fname3='test2_003.wav';
fname4='test2_004.wav';

nl=.1;
filname='ModulatedChirp_0.1SDutyTime_10.0kHzTo25.0kHz.wav';



[Ch1, Fs1] = audioread(fname1);
[Ch2, Fs2] = audioread(fname2);
[Ch3, Fs3] = audioread(fname3);
[Ch4, Fs4] = audioread(fname4);
[Fil,FilFs]=audioread(filname);

Ch1_noisy=Ch1+nl*randn(size(Ch1));
FS=Fs1;
time=1/FS*((1:length(Ch1))-1);


figure('position',[10,10,1200,800]);
subplot(2,2,1)
plot(time,Ch1,'linewidth',2)
hold on;
% plot(time,Ch2,'linewidth',2)
% plot(time,Ch3,'linewidth',2)
% plot(time,Ch4,'linewidth',2)
ax=gca;
ax.XLabel.String='Time (s)';
ax.YLabel.String='Amplitude (V)';
ax.FontSize=40;

subplot(2,2,2)
plot(time,Ch1_noisy,'linewidth',2)
hold on;
% plot(time,Ch2,'linewidth',2)
% plot(time,Ch3,'linewidth',2)
% plot(time,Ch4,'linewidth',2)
ax=gca;
ax.XLabel.String='Time (s)';
ax.YLabel.String='Amplitude (V)';
ax.FontSize=40;

[MF,lag]=xcorr(Ch1,Fil);
[MFnoisy,lagnoisy]=xcorr(Ch1_noisy,Fil);

subplot(2,2,3)
plot(MF)

subplot(2,2,4)
plot(MFnoisy)

figure('position',[0,0,1200,600])
plot(time,Ch1,'linewidth',2)
ax=gca;
ax.XLabel.String='Time (s)';
ax.YLabel.String='Amplitude (V)';
ax.FontSize=40;
exportgraphics(gcf,'/home/heathermoon/Desktop/Bootcamp/TankTest.png')

figure('position',[0,0,1200,600])
plot(time,Ch1_noisy,'linewidth',2)
ax=gca;
ax.XLabel.String='Time (s)';
ax.YLabel.String='Amplitude (V)';
ax.FontSize=40;
exportgraphics(gcf,'/home/heathermoon/Desktop/Bootcamp/TankTestwNoise.png')


figure('position',[0,0,1200,600])
plot(time,MF(lag>=0),'linewidth',2)
ax=gca;
ax.XLabel.String='Time (s)';
ax.YLabel.String='Amplitude (V)';
ax.FontSize=40;
exportgraphics(gcf,'/home/heathermoon/Desktop/Bootcamp/MatchFilterTankTest.png')

figure('position',[0,0,1200,600])
plot(time,MFnoisy(lagnoisy>=0),'linewidth',2)
ax=gca;
ax.XLabel.String='Time (s)';
ax.YLabel.String='Amplitude (V)';
ax.FontSize=40;
exportgraphics(gcf,'/home/heathermoon/Desktop/Bootcamp/MatchFilterTankTestwNoise.png')
