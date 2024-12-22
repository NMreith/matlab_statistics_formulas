close all; clear; clc,

x=zeros();
xd=[3312,14351,51496,77124,0,48212,33000,19048,19291,22,0,0,49712];

%minimum mumData 4;


%%Feed Data here @W and timeDiff, make sure there of even length.
%%if only one set of data is available copy it to both W and timeDiff
W=[3535,2914,2033,1743,1572,1412,1304,1209,1140,1077,1038,956,889,803,713,812,528,434,358,298,256,209,188,128,99,60,34,11,4,0];
timeDiff=[0.12,9/60,15/60,19/60,21/60,23/60,26/60,29/60,32/60,35/60,37/60,41/60,44/60,47/60,52/60,54/60,56/60,59/60,1,(1+(2/60)),(1+(3/60)),(1+(4/60)),(1+(5/60)),(1+(6/60)),(1+(8/60)),(1+(9/60)),(1+(10/60)),(1+(11/60)),(1+(12/60)),(1+(15/60))];

timeHr=[19,0];
n=size(timeDiff);
n_W=length(W);

n_T=length(timeDiff);
n_timeDiff=length(timeDiff);
xx=zeros();

var_W=zeros();
var_timeDiff=zeros();
R_h=zeros();

d_W=zeros();
d_W0=zeros();
d_timeDiff=zeros();
W_h=zeros();
timeDiff_h=zeros();

skwn_W=zeros();
skwn_W0=zeros();
krts_W=zeros();
krts_W0=zeros();
skwn_timeDiff=zeros();
krts_timeDiff=zeros();
timeDiffS=sort(timeDiff);

skwn_timeDiff0=zeros();
krts_timeDiff0=zeros();

std_W=zeros();
std_timeDiff=zeros();

SK_ai_W=zeros();
SK_af_W=zeros();
SK_b_W=zeros();

SK_ai_timeDiff=zeros();
SK_af_timeDiff=zeros();
SK_b_timeDiff=zeros();



m_W=((sum(W))./n_W);
m_timeDiff=((sum(timeDiff))./n_T);


%%%%%Mean Absolute Difference
mad_W0=(sum((W)-m_W));
mad_timeDiff=(sum((timeDiff)-m_timeDiff));


%%%%Variance
for i=1:n_W
    %
var_W0=var(W);
var_timeDiff0=var(timeDiff);
var_W(i)=(power((W(i)-m_W),2));
var_timeDiff(i)=(power((timeDiff(i)-m_timeDiff),2));
%var1=iarr./n;
%   iar=power(mad,2); 
end

%%%%Standard Deviation
for i=1:n_W
    std_W0=std(W);
%     std_timeDiff0=std(timeDiff);
    std_W(i)=sqrt(var_W(i));
    std_timeDiff0=sqrt(var_timeDiff0);
    std_timeDiff(i)=sqrt(var_timeDiff(i));
%    std1=stdd./n;
end

%%%%Median
%%Ungrouped
%
med_W=(n_W+1)/2;
med_timeDiff=(n_T+1)/2;
%k=class width
%i=((logd(n))/(logd(2)));

%%%%Mode
%%Ungrouped
%pick most repeated
%mod0=();
%
mod_i_W=(1);
mod_f_W=W(n_W);
mod_i_timeDiff=(1);
mod_f_timeDiff=timeDiff(n_T);

%%%%Harmonic Mean
%%Ungrouped Data
%%HM0=(n/(sum(1/xd)));
for i=1:n_W
    %
     W_h(i)=(1/W(i));
    HM0_W=(n_W/(sum(W_h)));
     timeDiff_h(i)=(1/timeDiff(i));
    HM0_timeDiff=(n_W/(sum(timeDiff_h)));
end

%%%%Measure of spread
%%Quatiles
%
Q1_W=((1/4).*(n_W+1));
Q1_W=W(round(Q1_W));
Q2_W=((2/4).*(n_W+1));
Q2_W=W(round(Q2_W));
Q3_W=((3/4).*(n_W+1));
Q3_W=W(round(Q3_W));

Q1_timeDiff=((1/4).*(n_T+1));
Q1_timeDiff=timeDiff(round(Q1_timeDiff));
Q2_timeDiff=((2/4).*(n_T+1));
Q2_timeDiff=timeDiff(round(Q2_timeDiff));
Q3_timeDiff=((3/4).*(n_T+1));
Q3_timeDiff=timeDiff(round(Q3_timeDiff));

%%Percentile
%
P1_W=((10/100).*(n_W+1));
P1_W=W(round(P1_W));
P5_W=((50/100).*(n_W+1));
P5_W=W(round(P5_W));
P9_W=((90/100).*(n_W));
P9_W=W(round(P9_W));

P1_timeDiff=((10/100).*(n_T+1));
P1_timeDiff=timeDiff(round(P1_timeDiff));
P5_timeDiff=((50/100).*(n_T+1));
P5_timeDiff=timeDiff(round(P5_timeDiff));
P9_timeDiff=((90/100).*(n_T));
P9_timeDiff=timeDiff(round(P9_timeDiff));




%%%%%%Measure of Skeewnesss and Kurtosis
%%  mean<median<mod -ve skewd
for i=1:1:n_W
    for j=1:1:n_W

SK_ai_W(i,j)=((m_W-mod_i_W)./std_W(j));
SK_af_W(i,j)=((m_W-mod_f_W)./std_W(j));
SK_b_W(i,j)=((3.*(m_W-med_W))./std_W(j));

SK_ai_timeDiff(i,j)=((m_W-mod_i_timeDiff)./std_timeDiff(j));
SK_af_timeDiff(i,j)=((m_timeDiff-mod_f_timeDiff)./std_timeDiff(j));
SK_b_timeDiff(i,j)=((3.*(m_timeDiff-med_timeDiff))./std_timeDiff(j));

    end 
end

%%%%%%Bowley's Co-efficient
%
Bc_W=((Q3_W-(2.*Q2_W)+Q1_W)./(Q3_W-Q1_W));
Bc_timeDiff=((Q3_timeDiff-(2.*Q2_timeDiff)+Q1_timeDiff)./(Q3_timeDiff-Q1_timeDiff));
%%10000000000000000000000000000000000000000000000000000000000000000000000000000
%%%%%%Kelly's Co-efficient
%
Kc_W=((P9_W-(2.*P5_W)+P1_W)./(P9_W-P1_W));
Kc_timeDiff=((P9_timeDiff-(2.*P5_timeDiff)+P1_timeDiff)./(P9_timeDiff-P1_timeDiff));


for i=1:n_W
    for j=0:1:n_W-1
         %
            d_W(i)=(W(i)-m_W);
            
                d_timeDiff(i)=(timeDiff(i)-m_timeDiff);
    
    for k=1:1:n(2)
        xx(k)=(power(10,k));
        x(k)=(k);
    end
    end
end
%%Skweaness and kurtosis
%
for i=1:1:n_W
    for j=1:1:n_W
        skwn_W(i,j)=((sum(power((d_W),3)))./(n_W.*(power(std_W(i),3))));
skwn_W0=((sum(power((d_W),3)))./(n_W.*(power(std_W0,3))));
krts_W(i)=((sum(power((d_W),4)))./(n_W.*(power(std_W(i),4))));
skwn_timeDiff(i)=((sum(power((d_timeDiff),3)))./(n_T.*(power(std_timeDiff(i),3))));
skwn_timeDiff0=((sum(power((d_timeDiff),3)))./(n_T.*(power(std_timeDiff0,3))));
krts_timeDiff0=((sum(power((d_timeDiff),4)))./(n_T.*(power(std_timeDiff0,4))));
krts_timeDiff(i)=((sum(power((d_timeDiff),4)))./(n_T.*(power(std_timeDiff(i),4))));

    end
end




%%%%Pearsons Product moment
Pr_W=(((((n_W.*(sum(W.*xx)))-(sum(W).*sum(xx))))./(sqrt((n_W.*(sum(power(W,2)))-(power((sum(W)),2))).*((n_W.*(sum(power(xx,2)))-(power((sum(W)),2))))))));
Pr_timeDiff=(((((n_T.*(sum(timeDiff.*xx)))-(sum(timeDiff).*sum(xx))))./(sqrt((n_T.*(sum(power(timeDiff,2)))-(power((sum(timeDiff)),2))).*((n_T.*(sum(power(xx,2)))-(power((sum(xx)),2))))))));
Pr0_W=(((((n_W.*(sum(W.*x)))-(sum(W).*sum(x))))./(sqrt((n_W.*(sum(power(W,2)))-(power((sum(W)),2))).*((n_W.*(sum(power(x,2)))-(power((sum(W)),2))))))));
Pr0_timeDiff=(((((n_T.*(sum(timeDiff.*x)))-(sum(timeDiff).*sum(x))))./(sqrt((n_T.*(sum(power(timeDiff,2)))-(power((sum(timeDiff)),2))).*((n_T.*(sum(power(xx,2)))-(power((sum(xx)),2))))))));
Pr_W_timediff=(((((n_W.*(sum(W.*timeDiff)))-(sum(W).*sum(timeDiff))))./(sqrt((n_W.*(sum(power(W,2)))-(power((sum(W)),2))).*((n_W.*(sum(power(timeDiff,2)))-(power((sum(W)),2))))))));

%%%%Regression Equation
%
%
b_W=(((((n_W.*(sum(W.*xx)))-(sum(W).*sum(xx))))./(sqrt((n_W.*(sum(power(W,2)))-(power((sum(W)),2)))))));
a_W=((sum(xx)-(b_W.*(sum(W)))));
b_timeDiff=(((((n_T.*(sum(timeDiff.*xx)))-(sum(timeDiff).*sum(xx))))./(sqrt((n_T.*(sum(power(timeDiff,2)))-(power((sum(timeDiff)),2)))))));
a_timeDiff=((sum(xx)-(b_timeDiff.*(sum(timeDiff)))));

b0_W=(((((n_W.*(sum(W.*x)))-(sum(W).*sum(x))))./(sqrt((n_W.*(sum(power(W,2)))-(power((sum(W)),2)))))));
a0_W=((sum(x)-(b_W.*(sum(W)))));
b0_timeDiff=(((((n_timeDiff.*(sum(timeDiff.*x)))-(sum(timeDiff).*sum(x))))./(sqrt((n_timeDiff.*(sum(power(timeDiff,2)))-(power((sum(timeDiff)),2)))))));
a0_timeDiff=((sum(x)-(b_timeDiff.*(sum(timeDiff)))));
%
b1_W=(((((n_W.*(sum(W.*timeDiff)))-(sum(W).*sum(timeDiff))))./(sqrt((n_W.*(sum(power(W,2)))-(power((sum(W)),2)))))));
a1_W=((sum(timeDiff)-(b_W.*(sum(W)))));


%%%%SpearMan Rank Correlation coefficient
%
rS_W=(1-((6.*(sum(power(d_W,2))))./(n_W.*(power(n,2)-1))));
rS_timeDiff=(1-((6.*(sum(power(d_timeDiff,2))))./(n_T.*(power(n,2)-1))));

figure
plot(d_timeDiff);
title('d_timeDiff');
figure
plot(d_W);
title('d_W');
figure
plot(krts_W);
title('krts_W');
figure
plot(krts_timeDiff);
title('krts_timeDiff');
figure
plot(SK_af_W);
title('SK_af_W');
figure
plot(SK_af_timeDiff);
title('SK_af_timeDiff');
figure
plot(SK_ai_W);
title('SK_ai_W');
figure
plot(SK_ai_timeDiff);
title('SK_ai_timeDiff');
figure
plot(SK_b_W);
title('SK_b_W');
figure
plot(SK_b_timeDiff);
title('SK_b_timeDiff');
figure
plot(skwn_W);
title('skwn_W');
figure
plot(skwn_timeDiff);
title('skwn_timeDiff');
figure
plot(std_W);
title('std_W');
figure
plot(std_timeDiff);
title('std_timeDiff');
figure
plot(timeDiff);
title('timeDiff');
figure
plot(timeDiffS);
title('timeDiffS');
figure
plot(var_W);
title('var_W');
figure
plot(var_timeDiff);
title('var_timeDiff');
figure
plot(W);
title('W');
figure
plot(W_h);
title('W_h');
figure
plot(timeDiff_h);
title('timeDiff_h');
figure
plot(x);
title('x');
figure
plot(xd);
title('xd');
figure
plot(xx);
title('xx');

%%en292-3045/2014
%%nelyi.maina.nm@gmail.com

