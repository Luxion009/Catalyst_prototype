clear Ine D_mat_mw D_mat_ss yy_dis zz_dis Imw_yy Imw_zz Imw_xx

			% %�������[�����g�v�Z
			% %	Ine(:,:)	=	Ixx	Ixy	Ixz		=	1	2	3
			% %					Iyx	Iyy	Iyz			4	5	6
			% %					Izx	Izy	Izz			7	8	9


			% disp("�嗃�ȊO�̊������[�����g����͂���")

			% Ine(1,1)=1914992314.76;
			% Ine(2,2)=9189816179.36;
			% Ine(3,3)=7276836877.38;
			% Ine(1,2)=-85542128.88;
			% Ine(2,1)=Ine(1,2);
			% Ine(1,3)=-3731721329.08;
			% Ine(3,1)=Ine(1,3);
			% Ine(2,3)=43873922.30;
			% Ine(3,2)=Ine(2,3);

			% Ine

			% %�������[�����g��(Kg/m^2)�ɕϊ�����
			% %10^3�������Ă��΂���
			% disp("unit converted");
			% Ine=Ine.*(10^(-9))


%�������[�����g��������
Ine(1:3,1:3)=0;


%�嗃�������[�����g�v�Z
%�������[�����g�Ɗ�����ς��v�Z����

%�񎟍\���̖��x������x�ɑ���
D_mat_ss(:,1)=Snd_rho.*S_mat(:,7);%�񎟍\�����ނ݂̂̐����x

D_mat_mw(:,1)=D_mat(:,1)./1000.+D_mat_ss(:,1);

%%%trapz�������Ă邩�킩��Ȃ����玩�O�Ōv�Z����
%���f�ԋ������o��
mat_dis(1,1)=S_mat(1,1);
for i=2:rows(S_mat(:,1))
	mat_dis(i,1)=S_mat(i,1)-S_mat(i-1,1);
endfor

%���f�ԋ��������ƂɊe��Ԃ̏d�ʂ����߂�
mat_mass(:,1)=mat_dis(:,1).*D_mat_mw(:,1);

wsmb_dis=1.01;	%input("�包���S���瓷�̃p�C�v���S�܂ł̋���(m)=");


%yy�������[�����g�悤�ɒ��S����̋������o��
yy_dis(:,1)=S_mat_swang(:,1).+wsmb_dis;

	%Imw_yy=trapz(S_mat(:,1),D_mat_mw(:,1).*(yy_dis(:,1)).^2);
Imw_yy=sum(yy_dis(:,1).^2.*mat_mass(:,1));
Imw_yy=Imw_yy.*2

			% yz�悤��y*z������
			% 	yz_dis(:,1)=S_mat(:,1).*(S_mat_swang(:,1).+wsmb_dis);

			% 	Imw_yz=trapz(yz_dis(:,1),D_mat_mw(:,1).*yz_dis(:,1));
			% 	Imw_yz=Imw_yz.*2
			% �~�X�A�Ώ̂����犵����ς̓[��


%zz
zz_dis(:,1)=S_mat(:,1);
	%Imw_zz=trapz(S_mat(:,1),D_mat_mw(:,1).*(zz_dis(:,1)).^2);
Imw_zz=sum(zz_dis(:,1).^2.*mat_mass(:,1));		
Imw_zz=Imw_zz.*2

%xx
xx_dis(:,1)=sqrt(S_mat(:,1).^2.+(S_mat_swang(:,1).+wsmb_dis).^2);
Imw_xx=sum(xx_dis(:,1).^2.*mat_mass(:,1));
Imw_xx=Imw_xx.*2

			% %�p�C���b�g�������[�����g�v�Z
			% pw_dis=0.121;	%input("�d�S����p�C���b�g�d�S�܂ł̋���(m)=");
			% p_mass=66;	%input("�p�C���b�g�̑̏d(kg)=");

			% Ipw_yy=p_mass*pw_dis^2;
			% Ipw_xx=Ipw_yy

			% %�y���������[�����g���v�Z
			% r_mass=1.2;
			% rw_dis=1.75;

			% Irw_yy=rw_dis^2*r_mass;
			% Irw_zz=Irw_yy


%�e�������[�����g�����Z����
%�y��
Ine_per=calc_ine(1.4,1.797,0,0.805)
%��쓮��
Ine_upd=calc_ine(0.65,1.399,0,0.799)
%���쓮��
Ine_lwd=calc_ine(1.16,0.919,0,-0.351)
%�O�^�C��
Ine_ft=calc_ine(0.3,0.866,0,-0.757)
%���^�C��
Ine_bt=calc_ine(0.3,-0.583,0,-0.757)
%�t���[��
Ine_fl=calc_ine(3.95,-0.92,0,0.381)
%�֎q�ƃp�C���b�g�ƃt�F�A�����O
Ine_pa=calc_ine(67.9,0.118,0,-0.387)
%�e�[���p�C�v
Ine_tl=calc_ine(1.4,-2.955,0,0.798)
%��������
Ine_ht=calc_ine(0.75,-4.513,0,0.909)
%��������
Ine_vt=calc_ine(0.75,-5.079,0,0.849)

%�S�������[�����g�ɉ��Z����

Ine(1,1)=Ine(1,1).+Imw_xx;
Ine(2,2)=Ine(2,2).+Imw_yy;
Ine(3,3)=Ine(3,3).+Imw_zz;

Ine=Ine+Ine_per+Ine_upd+Ine_lwd+Ine_ft+Ine_bt+Ine_fl+Ine_pa+Ine_tl+Ine_ht+Ine_vt