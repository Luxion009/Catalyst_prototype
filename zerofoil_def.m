%�������[���̗��^��Cl*C�𓾂�֐�
function clc=zerofoil_def(cl_c,aoa,vel,adata)
nu=1.604*10^(-5);
c_list=0.2:0.1:0.8;
ang_list=-3:1:7;
re_list=c_list.*vel/nu;

disp("generating zerofoil");

for rec=1:length(c_list)
	printf("%d,",rec);
	% for angc=1:11
	% 	angc
	% 	% angc_cl(angc,1)=Re_spline(re_list(rec),6,1,ang_list(angc),2,adata);
	% 	angc_cl(angc,1)=xfoil(adata,ang_list(angc),re_list(rec),0).CL;
		
	% endfor
	
	% a_cl(1,rec)=spline(ang_list,angc_cl(:,1),aoa);
	a_cl(1,rec)=xfoil(adata,aoa,re_list(rec),0).CL;
	%����RE���ɂ����邠��AOA��CL
	
endfor

disp("")

%CL�����Ƃ�C��CL���v�Z����
for i=1:length(c_list)
	c_cl_o(1,i)=a_cl(1,i)*c_list(1,i);
endfor

%���͂���c*cl�ɂ���C��spline�ł��Ƃ߂�
% length(c_cl_o)
% length(c_list)

clc=interp1(c_cl_o,c_list,cl_c,"spline");

endfunction
