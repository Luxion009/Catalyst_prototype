%垂直尾翼定義

lvd=5.089;

Vv=input("垂直尾翼容積=");
Sv=Vv*(span+wl_span*2)*mw_s/lvd
rs_v=input("リブ間=");
rnu_v=input("上部リブ間数=");
rnm_v=input("中央リブ間数=");
rnl_v=input("下部リブ間数=");

chord_ev=input("端部コード長(m)=");

chord_cv=(Sv-0.5*rs_v*(rnu_v+rnl_v)*chord_ev)/((0.5*rs_v*(rnu_v+rnl_v+2*rnm_v)))%中央部コード長を求める

span_v=rs_v*(rnu_v+rnm_v+rnl_v)
vt_aspect=span_v^2/Sv
Fv=Sv*lvd^2/mw_s/((span+wl_span*2)^2)

%垂直揚力傾斜を近似で求める
%dClv=7.448/(1+7.448/(pi*(vt_aspect)));%0010
%dClv=6.417/(1+6.417/(pi*(vt_aspect)));%0009
dClv=6.245/(1+6.245/(pi*(vt_aspect)));%sd8020
%dClv=3.7581;

Clbetafin=0;%胴体の上下に出る形だから無視できるほど小さいと過程
Cybeta=-Sv/mw_s*dClv;
Cyr=Sv/mw_s*dClv*2*lvd/(span+wl_span*2);
Cnbetafin=Vv*dClv;
Cndr=-Vv*dClv;
Cnrfin=-Vv*dClv*2*lvd/(span+wl_span*2);

x_dcoef=a_density*d_vel^2*mw_s*(span+wl_span*2)/2/Ine(1,1);
y_dcoef=a_density*d_vel^2*mw_s*(span+wl_span*2)/2/Ine(2,2);
z_dcoef=a_density*d_vel^2*mw_s*(span+wl_span*2)/2/Ine(3,3);

Cybeta;
Clbeta;
Cnbetaall=Cnbeta+Cnbetafin;
Clp;
Cnp;
Clr;
Cnr;

Lbeta=x_dcoef*Clbeta
Lp=x_dcoef*Clp/2*(span+wl_span*2)/d_vel
Lr=x_dcoef*Clr/2*(span+wl_span*2)/d_vel

Ybeta=a_density*d_vel*mw_s/2/L_lst*Cybeta;
Yr=a_density*d_vel*mw_s*(span+wl_span*2)/2/L_lst*Cyr;

Nbeta=z_dcoef*Cnbetaall
Np=z_dcoef*Cnp/2*(span+wl_span*2)/d_vel
Nr=z_dcoef*(Cnr+Cnrfin)/2*(span+wl_span*2)/d_vel
Ndr=z_dcoef*Cndr;

printf('スパイラルモード\nLbeta*Nr-Nbeta*Lr:%f\n時定数:%f\n',Lbeta*Nr-Nbeta*Lr,(-Lp+Lbeta/Nbeta*(Np-9.8/d_vel))/(9.8/d_vel*(Lbeta/Nbeta*Nr-Lr)));