%% Calculate m and SD from mu and sigma

clear
clc
GMMTab=readtable("C:\Users\lasss\Documents\Research\Brewer Lab work\All Presentations\Sam\Spike Field Coherence\byAxonGMMFits 260707b.xlsx");

mu1=GMMTab.mu2;
sigma1=GMMTab.sigma2;

[m1,v1]=lognstat(mu1,sigma1);
SD1=sqrt(sigma1);

GMMTab.m1=m1;
GMMTab.SD1=SD1;

myTab=GMMTab(:,["fi","Subregion","Electrode","m1","SD1"]);

%% Dunnett

mu1ECDG=GMMTab.mu1(GMMTab.Subregion=="EC-DG");
mu2ECDG=GMMTab.mu2([8:12]);
gECDG=[repmat("ECDG",1,sum(GMMTab.Subregion=="EC-DG")),"3F4","1F2","6F4","5F2","1F4"];
figure
[~,~,stats]=anovan([mu1ECDG;mu2ECDG],{gECDG},"display","off");
[c,~]=multcompare(stats,"CriticalValueType","dunnett");

mu1DGCA3=GMMTab.mu1(GMMTab.Subregion=="DG-CA3");
mu2DGCA3=GMMTab.mu2([29:38]);
gDGCA3=[repmat("DGCA3",1,sum(GMMTab.Subregion=="DG-CA3")),"3J6","2M6","3K6","4L6","1K6","6J6","4K6","5M6","5L6","4M6"];
figure
[~,~,stats]=anovan([mu1DGCA3;mu2DGCA3],{gDGCA3},"display","off");
[c,~]=multcompare(stats,"CriticalValueType","dunnett");

mu1CA3CA1=GMMTab.mu1(GMMTab.Subregion=="CA3-CA1");
mu2CA3CA1=GMMTab.mu2([47:52]);
gCA3CA1=[repmat("CA3CA1",1,sum(GMMTab.Subregion=="CA3-CA1")),"1G11","3G10","3G12","5G12","5G11","6G12"];
figure
[~,~,stats]=anovan([mu1CA3CA1;mu2CA3CA1],{gCA3CA1},"display","off");
[c,~]=multcompare(stats,"CriticalValueType","dunnett");

mu1CA1EC=GMMTab.mu1(GMMTab.Subregion=="CA1-EC");
mu2CA1EC=GMMTab.mu2([61:64]);
gCA1EC=[repmat("CA1EC",1,sum(GMMTab.Subregion=="CA1-EC")),"4C7","1D7","4A7","6C7"];
figure
[~,~,stats]=anovan([mu1CA1EC;mu2CA1EC],{gCA1EC},"display","off");
[c,~]=multcompare(stats,"CriticalValueType","dunnett");

mu1ECCA3=GMMTab.mu1(GMMTab.Subregion=="EC-CA3");
mu2ECCA3=GMMTab.mu2([19:23]);
gECCA3=[repmat("ECCA3",1,sum(GMMTab.Subregion=="EC-CA3")),"4E7","4F6","6F6","4F5","6G5"];
figure
[~,~,stats]=anovan([mu1ECCA3;mu2ECCA3],{gECCA3},"display","off");
[c,~]=multcompare(stats,"CriticalValueType","dunnett");