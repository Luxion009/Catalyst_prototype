%X�����Ȃ��v�Z
for n=1:div_c %��������
	%-----���W�A�ϑw�f�[�^�A�}���h�����f�[�^�A���ԍ��A�����𑗂��Ċe���W�_�ł̒f�ʓ񎟃��[�����g���v�Z����֐�(kgf.m^2)
	%disp("generating S_I");
	S_I(dp_num*(n-1)+1:dp_num*n,1:2)=calc_i(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,E,C_T)./1000000;
	%-----���W�A�ϑw�f�[�^�A�}���h�����f�[�^�A���ԍ��A�����𑗂��Ċe���W�_�ł̋Ȃ��������v�Z����֐�(kgf.m^2)
	%disp("generating S_H");
	S_H(dp_num*(n-1)+1:dp_num*n,1:2)=calc_hardness(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,E,C_T)./1000000;
	%-----�e���W�ł̐����x���v�Z����֐�(kg/m)
	%disp("generating D_mat");
	if n!=div_c
		D_mat(dp_num*(n-1)+1:dp_num*n,1)=calc_condens(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},ply_data{n+1},mand_data,n,C_D,C_T).*1000;
	else
		D_mat(dp_num*(n-1)+1:dp_num*n,1)=calc_ldens(M_data(dp_num*(n-1)+1:dp_num*n,1).*1000,ply_data{n},mand_data,n,C_D,C_T).*1000;

	endif
endfor	

%�ޗ��v�Z�ɕK�v�ȃf�[�^��ʕϐ��Ɋi�[
disp("Data transfer");
S_mat(:,1)=M_data(:,1);
S_mat(:,6)=D_mat(:,1);%���W�ɂ���������x
S_mat(:,7)=M_data(:,2);%�R�[�h��
S_mat(:,2)=M_data(:,3);%�Ǐ��g�͌W��
S_mat(:,3)=M_data(:,4);%�Ǐ��R�͌W��
S_mat(:,4)=M_data(:,5)./180.*pi;%�U���p�x�irad�j
S_mat(:,5)=M_data(:,6);%�������x
ID_swang(:,1)=interp1(nang_line(1:250,1),line_e(1:250,2),S_mat(:,1),'spline','extrap');%���񂾎���Z�������W

disp("genarating d_D");
%-----�Ǐ��R�͕��z(N/m)
S_mat_dD(:,1)=(S_mat(:,3).*cos(S_mat(:,4))+S_mat(:,2).*sin(S_mat(:,4))).*0.5.*S_mat(:,7).*a_density.*(S_mat(:,5).^2);

%disp("generating Q");
for i=1:div_c*dp_num-1
	S_mat_dDQ(i,1)=(S_mat_dD(i,1)+S_mat_dD(i+1,1))/2*(S_mat(i+1,1)-S_mat(i,1));
end
S_mat_dDQ(div_c*dp_num,1)=interp1(S_mat(1:div_c*dp_num-1,1),S_mat_dDQ(:,1),S_mat(div_c*dp_num,1),'linear','extrap');

S_mat_DQ(div_c*dp_num+1,1)=0;
for i=div_c*dp_num:-1:1
	S_mat_DQ(i,1)=S_mat_DQ(i+1,1)+S_mat_dDQ(i,1);
end
%-----����f��(kgf)
S_mat_DQ=S_mat_DQ./9.8;

%disp("generating M");
for i=1:div_c*dp_num-1
	S_mat_dDM(i,1)=(S_mat_DQ(i,1)+S_mat_DQ(i+1,1))/2*(S_mat(i+1,1)-S_mat(i,1))*1000;
end
%-----���[�����g(kgf.m)
S_matDM(div_c*dp_num,1)=0;
for i=div_c*dp_num-1:-1:1
	S_matDM(i,1)=S_matDM(i+1,1)+S_mat_dDM(i,1);
end

%disp("generating Xx");
for n=1:div_c
	Xx(dp_num*(n-1)+1:dp_num*n,1)=S_I(dp_num*(n-1)+1:dp_num*n,1)./ply_data{n}(rows(ply_data{n}),1).*1000;%�f�ʌW��
end

for n=1:div_c
	S_max_X(1:dp_num*n,1)=S_matDM(1:dp_num*n,1)./Xx(1:dp_num*n,1).*9.8./1000;%�ő�Ȃ�����

end

for i=1:div_c*dp_num-1
	S_mat_dDtheta(i,1)=(S_matDM(i,1)/S_H(i,1)+S_matDM(i+1,1)/S_H(i+1,1))/2*(S_mat(i+1,1)-S_mat(i,1));%�A�������̋ȗ��̕��ςɂ��̊Ԃ̋�����������
	%�ڐ��̊p�x�ω�
end
S_mat_Dtheta(1,1)=0;%�ڐ��Ɛ������̊p�x
for i=2:div_c*dp_num-1
	S_mat_Dtheta(i,1)=S_mat_Dtheta(i-1,1)+S_mat_dDtheta(i,1);
end
S_mat_Dtheta(div_c*dp_num,1)=interp1(S_mat(1:div_c*dp_num-1,1),S_mat_Dtheta(:,1),S_mat(div_c*dp_num,1),'linear','extrap');

for i=1:div_c*dp_num-1
	S_mat_dDswang(i,1)=(S_mat_Dtheta(i,1)+S_mat_Dtheta(i+1,1))./2.*(S_mat(i+1,1)-S_mat(i,1));%sin x=x�Ƌߎ��ł���̂�X�����ω����Z�o�ł���
end
S_mat_Dswang(1,1)=S_mat_dDswang(1,1);
for i=2:div_c*dp_num-1
	S_mat_Dswang(i,1)=S_mat_Dswang(i-1,1)+S_mat_dDswang(i,1);%X�����ʒu
end
S_mat_Dswang(div_c*dp_num,1)=interp1(S_mat(1:div_c*dp_num-1,1),S_mat_Dswang(:,1),S_mat(div_c*dp_num,1),'linear','extrap');
S_mat_Dswang(:,1)=S_mat_Dswang.*10^(-3);