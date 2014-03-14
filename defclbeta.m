%�����ׂ�̈�����W�����݂���
%�㔼�p����

clear dw_ang_cp dw_ang ds_ang

slipang=3 %������p(deg)

%�Ǐ��㔼�p�����߂�(500�_�(rad))
dw_ang(1,1)=0;
for i=2:rows(S_mat(:,1))-1
	dw_ang(i,1)=atan((S_mat_swang(i+1,1)-S_mat_swang(i-1,1))/(S_mat(i+1,1)-S_mat(i-1,1)));
endfor
dw_ang(rows(S_mat(:,1)),1)=atan((S_mat_swang(rows(S_mat(:,1)),1)-S_mat_swang(rows(S_mat(:,1))-1,1))/(S_mat(rows(S_mat(:,1)),1)-S_mat(rows(S_mat(:,1))-1,1)));

%�Ǐ��㔼�p��CP�ɕ⊮����(rad)
dw_ang_cp(:,1)=interp1(S_mat(:,1),dw_ang(:,1),nang_cp(1:250,1),'spline','extrap');
dw_ang_cp(251:270,1)=dw_ang_cp(250,1)-(wl_ang*pi/180);%500�_�f�[�^�̒��ɃE�C���O���b�g�͊܂܂�Ȃ�����


%�����}�p(deg)
ds_ang=atan((d_vel.*sin(slipang*pi/180).*sin(dw_ang_cp(:,1)))./d_vel).*180./pi;

%wing_analysis

%%%���񎞓���
disp("calculating innerwing")

%��`�����f�[�^���R���g���[���|�C���g�ɕۊǂ���

lst_fdata(:,1)=nang_cp(:,1);
lst_fdata(:,2)=interp1(nang_line(:,1),f_data(:,2),nang_cp(:,1),"extrap");
lst_fdata(:,3)=interp1(nang_line(:,1),f_data(:,3),nang_cp(:,1),"extrap");

for i=1:50
	lst_ydata(:,i)=y_foil(:,1);
endfor

for i=2:w_div_c-1
	for j=1:dw_divn
		lst_ydata(:,(i-1)*dw_divn+j)=y_foil(:,i-1)*lst_fdata(i,3)+y_foil(:,i)*(1-lst_fdata(i,3));
	endfor
endfor

for j=1:dw_divn
	lst_ydata(:,4*dw_divn+j)=y_foil(:,4);
endfor

for j=we_div+1:we_div+wld_c
	lst_ydata(:,j)=y_foil(:,4);
endfor

lst_gamma=Hex_gamma;


if(c!=1)
	disp("regenerting lst_gamma");
	lst_gamma(:,1)=lst_data(:,1).*lst_fdata(:,2).*lst_comp_vel(:,1)./2;
endif

%���������U�����x
disp("recalculating Vn_i.");
for i=1:(wld_c+we_div)
	for j=1:(wld_c+we_div)
		lst_Vn_i_d(j,i)=Q_ij(i,j)*lst_gamma(j,1)./2;
	endfor
	lst_Vn_i(i,1)=sum(lst_Vn_i_d(:,i));
	%printf("%d/100\n",i);
endfor

%�U���p�x
disp("calculating i_ang.");
for i=1:(wld_c+we_div)
	lst_i_ang(i,1)=atan(lst_Vn_i(i,1)/d_vel)*180/pi;
endfor

%�������x�v�Z
disp("calculating comp_vel");
for i=1:(wld_c+we_div)
	lst_comp_vel(i,1)=sqrt(d_vel^2+lst_Vn_i(i,1)^2);
endfor

disp("generating lst_data");


errc=1;
for n=1:we_div+wld_c
	lst_re(n,1)=lst_fdata(n,2)*lst_comp_vel(n,1)/nu;
	[bufx foil err]=xfoil([x_foil lst_ydata(:,n)],f_mtang(n,1)-lst_i_ang(n,1)+ds_ang(n,1),lst_re(n,1),0);
	if(err==0)
		printf("%d,",n);
		lst_data(n,:)=[bufx.CL bufx.CD bufx.Cm];
	else
		printf("e,");
		errlist(errc,1)=n;
		errc++;
	endif
end

if errc!=1
	for n=1:errc-1
		ncc=errlist(n,1)
		lst_data(ncc,:)=[interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,1),nang_cp(ncc,1)) interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,2),nang_cp(ncc,1)) interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,3),nang_cp(ncc,1))];
	endfor
endif
buff=input("")

disp("regenerting lst_gamma");
lst_gamma(:,1)=lst_data(:,1).*lst_fdata(:,2).*lst_comp_vel(:,1)./2;

%���������U�����x
disp("recalculating Vn_i.");
for i=1:(wld_c+we_div)
	for j=1:(wld_c+we_div)
		lst_Vn_i_d(j,i)=Q_ij(i,j)*lst_gamma(j,1)./2;
	endfor
	lst_Vn_i(i,1)=sum(lst_Vn_i_d(:,i));
	%printf("%d/100\n",i);
endfor

%�U���p�x
disp("calculating i_ang.");
for i=1:(wld_c+we_div)
	lst_i_ang(i,1)=atan(lst_Vn_i(i,1)/d_vel)*180/pi;
endfor

%�������x�v�Z
disp("calculating comp_vel");
for i=1:(wld_c+we_div)
	lst_comp_vel(i,1)=sqrt(d_vel^2+lst_Vn_i(i,1)^2);
endfor

disp("");
disp("generating L_lst_d");

for i=1:we_div
	L_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,1).*cos(lst_i_ang(i,1)./180.*pi)-lst_data(i,2).*sin(lst_i_ang(i,1)./180.*pi)).*cos(line_e_ang(i,1));

endfor

for i=we_div:we_div+wld_c
	L_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,1).*cos(lst_i_ang(i,1)./180.*pi)-lst_data(i,2).*sin(lst_i_ang(i,1)./180.*pi)).*cos(line_e_ang(i,1));

endfor

%�ߎ���
for i=1:5
	L_lst_d(i*dw_divn-1:i*dw_divn,1)=L_lst_d(i*dw_divn-2,1);
	L_lst_d(i*dw_divn+1:i*dw_divn+2,1)=L_lst_d(i*dw_divn+3,1);			
endfor

%�������[��������
dM_ire(:,1)=L_lst_d(:,1);

for i=1:we_div
	Dp_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,2).*cos(lst_i_ang(i,1)./180.*pi)+lst_data(i,1).*sin(lst_i_ang(i,1)./180.*pi));
	
endfor

for i=we_div:we_div+wld_c
	Dp_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,2).*cos(lst_i_ang(i,1)./180.*pi)+lst_data(i,1).*sin(lst_i_ang(i,1)./180.*pi));
	
endfor

%�ߎ���
for i=1:5
	Dp_lst_d(i*dw_divn-1:i*dw_divn,1)=Dp_lst_d(i*dw_divn-2,1);
	Dp_lst_d(i*dw_divn+1:i*dw_divn+2,1)=Dp_lst_d(i*dw_divn+3,1);			
endfor

%�������[�_���s���O
dM_ire(:,2)=Dp_lst_d(:,1);

%%%���񎞊O��
disp("calculating outerwing")

%��`�����f�[�^���R���g���[���|�C���g�ɕۊǂ���

lst_fdata(:,1)=nang_cp(:,1);
lst_fdata(:,2)=interp1(nang_line(:,1),f_data(:,2),nang_cp(:,1),"extrap");
lst_fdata(:,3)=interp1(nang_line(:,1),f_data(:,3),nang_cp(:,1),"extrap");

for i=1:50
	lst_ydata(:,i)=y_foil(:,1);
endfor

for i=2:w_div_c-1
	for j=1:dw_divn
		lst_ydata(:,(i-1)*dw_divn+j)=y_foil(:,i-1)*lst_fdata(i,3)+y_foil(:,i)*(1-lst_fdata(i,3));
	endfor
endfor

for j=1:dw_divn
	lst_ydata(:,4*dw_divn+j)=y_foil(:,4);
endfor

for j=we_div+1:we_div+wld_c
	lst_ydata(:,j)=y_foil(:,4);
endfor

lst_gamma=Hex_gamma;


if(c!=1)
	disp("regenerting lst_gamma");
	lst_gamma(:,1)=lst_data(:,1).*lst_fdata(:,2).*lst_comp_vel(:,1)./2;
endif

%���������U�����x
disp("recalculating Vn_i.");
for i=1:(wld_c+we_div)
	for j=1:(wld_c+we_div)
		lst_Vn_i_d(j,i)=Q_ij(i,j)*lst_gamma(j,1)./2;
	endfor
	lst_Vn_i(i,1)=sum(lst_Vn_i_d(:,i));
	%printf("%d/100\n",i);
endfor

%�U���p�x
disp("calculating i_ang.");
for i=1:(wld_c+we_div)
	lst_i_ang(i,1)=atan(lst_Vn_i(i,1)/d_vel)*180/pi;
endfor

%�������x�v�Z
disp("calculating comp_vel");
for i=1:(wld_c+we_div)
	lst_comp_vel(i,1)=sqrt(d_vel^2+lst_Vn_i(i,1)^2);
endfor

disp("generating lst_data");


errc=1;
for n=1:we_div+wld_c
	lst_re(n,1)=lst_fdata(n,2)*lst_comp_vel(n,1)/nu;
	[bufx foil err]=xfoil([x_foil lst_ydata(:,n)],f_mtang(n,1)-lst_i_ang(n,1)-ds_ang(n,1),lst_re(n,1),0);
	if(err==0)
		printf("%d,",n);
		lst_data(n,:)=[bufx.CL bufx.CD bufx.Cm];
	else
		printf("e,");
		errlist(errc,1)=n;
		errc++;
	endif
end

if errc!=1
	for n=1:errc-1
		ncc=errlist(n,1)
		lst_data(ncc,:)=[interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,1),nang_cp(ncc,1)) interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,2),nang_cp(ncc,1)) interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,3),nang_cp(ncc,1))];
	endfor
endif
buff=input("")

disp("regenerting lst_gamma");
lst_gamma(:,1)=lst_data(:,1).*lst_fdata(:,2).*lst_comp_vel(:,1)./2;

%���������U�����x
disp("recalculating Vn_i.");
for i=1:(wld_c+we_div)
	for j=1:(wld_c+we_div)
		lst_Vn_i_d(j,i)=Q_ij(i,j)*lst_gamma(j,1)./2;
	endfor
	lst_Vn_i(i,1)=sum(lst_Vn_i_d(:,i));
	%printf("%d/100\n",i);
endfor

%�U���p�x
disp("calculating i_ang.");
for i=1:(wld_c+we_div)
	lst_i_ang(i,1)=atan(lst_Vn_i(i,1)/d_vel)*180/pi;
endfor

%�������x�v�Z
disp("calculating comp_vel");
for i=1:(wld_c+we_div)
	lst_comp_vel(i,1)=sqrt(d_vel^2+lst_Vn_i(i,1)^2);
endfor


disp("");
disp("generating L_lst_d");

for i=1:we_div
	L_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,1).*cos(lst_i_ang(i,1)./180.*pi)-lst_data(i,2).*sin(lst_i_ang(i,1)./180.*pi)).*cos(line_e_ang(i,1));

endfor

for i=we_div:we_div+wld_c
	L_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,1).*cos(lst_i_ang(i,1)./180.*pi)-lst_data(i,2).*sin(lst_i_ang(i,1)./180.*pi)).*cos(line_e_ang(i,1));

endfor

%�ߎ���
for i=1:5
	L_lst_d(i*dw_divn-1:i*dw_divn,1)=L_lst_d(i*dw_divn-2,1);
	L_lst_d(i*dw_divn+1:i*dw_divn+2,1)=L_lst_d(i*dw_divn+3,1);			
endfor

%�O�����[��������
dM_ore(:,1)=L_lst_d(:,1);

for i=1:we_div
	Dp_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,2).*cos(lst_i_ang(i,1)./180.*pi)+lst_data(i,1).*sin(lst_i_ang(i,1)./180.*pi));
	
endfor

for i=we_div:we_div+wld_c
	Dp_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,2).*cos(lst_i_ang(i,1)./180.*pi)+lst_data(i,1).*sin(lst_i_ang(i,1)./180.*pi));
	
endfor

%�ߎ���
for i=1:5
	Dp_lst_d(i*dw_divn-1:i*dw_divn,1)=Dp_lst_d(i*dw_divn-2,1);
	Dp_lst_d(i*dw_divn+1:i*dw_divn+2,1)=Dp_lst_d(i*dw_divn+3,1);			
endfor

%�O�����[�_���s���O
dM_ore(:,2)=Dp_lst_d(:,1);

%�͂����[�����g��
M_ire(1,1)=trapz(lst_fdata(:,1),dM_ire(:,1).*lst_fdata(:,1));
M_ire(2,1)=trapz(lst_fdata(:,1),dM_ire(:,2).*lst_fdata(:,1));

M_ore(1,1)=trapz(lst_fdata(:,1),dM_ore(:,1).*lst_fdata(:,1));
M_ore(2,1)=trapz(lst_fdata(:,1),dM_ore(:,2).*lst_fdata(:,1));

%�������[�����g
re_m=M_ire(1,1)-M_ore(1,1);
%���[�_���s���O
re_y=M_ire(2,1)-M_ore(2,1);

Clbeta=-re_m/Kp/span/slipang*180/pi-0.0004
Cnbeta=re_y/Kp/span/slipang*180/pi