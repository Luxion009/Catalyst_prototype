%����������`

lvd=5.089;

Vv=input("���������e��=");
Sv=Vv*(span+wl_span*2)*mw_s/lvd
rs_v=input("���u��=");
rnu_v=input("�㕔���u�Ԑ�=");
rnm_v=input("�������u�Ԑ�=");
rnl_v=input("�������u�Ԑ�=");

chord_ev=input("�[���R�[�h��(m)=");

chord_cv=(Sv-0.5*rs_v*(rnu_v+rnl_v)*chord_ev)/((0.5*rs_v*(rnu_v+rnl_v+2*rnm_v)))%�������R�[�h�������߂�

span_v=rs_v*(rnu_v+rnm_v+rnl_v)
vt_aspect=span_v^2/Sv
Fv=Sv*lvd^2/mw_s/((span+wl_span*2)^2)

%�����g�͌X�΂��ߎ��ŋ��߂�
%dClv=7.448/(1+7.448/(pi*(vt_aspect)));%0010
%dClv=6.417/(1+6.417/(pi*(vt_aspect)));%0009
dClv=6.245/(1+6.245/(pi*(vt_aspect)));%sd8020
%dClv=3.7581;

Clbetafin=0;%���̂̏㉺�ɏo��`�����疳���ł���قǏ������Ɖߒ�
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

printf('�X�p�C�������[�h\nLbeta*Nr-Nbeta*Lr:%f\n���萔:%f\n',Lbeta*Nr-Nbeta*Lr,(-Lp+Lbeta/Nbeta*(Np-9.8/d_vel))/(9.8/d_vel*(Lbeta/Nbeta*Nr-Lr)));