%�e�������ɂ�����f�[�^�]���_���ꗃ�P�O�O�_�ɑ������
dp_num=100;
div_c=w_div_c;

%���W�_���`
disp("generating line element");
cw=1;
cw
for i=1:dp_num
	M_data(i,1)=divw_span(1,1)./dp_num.*i;
endfor

for cw=2:w_div_c
	cw

	for i=1:dp_num
		M_data((cw-1)*dp_num+i,1)=f_data(dw_divn*(cw-1),1)+divw_span(cw,1)./dp_num.*i;
	endfor	
endfor

%�R���g���[���|�C���g�̉�̓f�[�^����f�[�ōĒ�`����
lst_data_ed(:,1)=interp1(nang_cp(:,1),lst_data(:,1),nang_line(:,1),"extrap");
lst_data_ed(:,2)=interp1(nang_cp(:,1),lst_data(:,2),nang_line(:,1),"extrap");
lst_i_ang_ed(:,1)=interp1(nang_cp(:,1),lst_i_ang(:,1),nang_line(:,1),"extrap");
lst_comp_vel_ed(:,1)=interp1(nang_cp(:,1),lst_comp_vel(:,1),nang_line(:,1),"extrap");
disp("Foildata was converted.");


%���W�X�����ƂɃf�[�^��⊮
for i=1:w_div_c*dp_num
	M_data(i,2)=interp1(f_data(:,1),f_data(:,2),M_data(i,1),'spline','extrap');
	M_data(i,3)=interp1(f_data(:,1),lst_data_ed(:,1),M_data(i,1),'spline','extrap');
	M_data(i,4)=interp1(f_data(:,1),lst_data_ed(:,2),M_data(i,1),'spline','extrap');
	M_data(i,5)=interp1(f_data(:,1),lst_i_ang_ed(:,1),M_data(i,1),'spline','extrap');
	M_data(i,6)=interp1(f_data(:,1),lst_comp_vel_ed(:,1),M_data(i,1),'spline','extrap');
endfor

disp("Foildata was extended.");