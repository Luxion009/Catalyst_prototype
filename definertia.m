clear Ine D_mat_mw D_mat_ss yy_dis zz_dis Imw_yy Imw_zz Imw_xx

			% %慣性モーメント計算
			% %	Ine(:,:)	=	Ixx	Ixy	Ixz		=	1	2	3
			% %					Iyx	Iyy	Iyz			4	5	6
			% %					Izx	Izy	Izz			7	8	9


			% disp("主翼以外の慣性モーメントを入力する")

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

			% %慣性モーメントを(Kg/m^2)に変換する
			% %10^3をかけてやればいい
			% disp("unit converted");
			% Ine=Ine.*(10^(-9))


%慣性モーメントを初期化
Ine(1:3,1:3)=0;


%主翼慣性モーメント計算
%慣性モーメントと慣性乗積を計算する

%二次構造の密度を線密度に足す
D_mat_ss(:,1)=Snd_rho.*S_mat(:,7);%二次構造部材のみの線密度

D_mat_mw(:,1)=D_mat(:,1)./1000.+D_mat_ss(:,1);

%%%trapzがあってるかわからないから自前で計算する
%線素間距離を出す
mat_dis(1,1)=S_mat(1,1);
for i=2:rows(S_mat(:,1))
	mat_dis(i,1)=S_mat(i,1)-S_mat(i-1,1);
endfor

%線素間距離をもとに各区間の重量を求める
mat_mass(:,1)=mat_dis(:,1).*D_mat_mw(:,1);

wsmb_dis=1.01;	%input("主桁中心から胴体パイプ中心までの距離(m)=");


%yy慣性モーメントように中心からの距離を出す
yy_dis(:,1)=S_mat_swang(:,1).+wsmb_dis;

	%Imw_yy=trapz(S_mat(:,1),D_mat_mw(:,1).*(yy_dis(:,1)).^2);
Imw_yy=sum(yy_dis(:,1).^2.*mat_mass(:,1));
Imw_yy=Imw_yy.*2

			% yzようにy*zをだす
			% 	yz_dis(:,1)=S_mat(:,1).*(S_mat_swang(:,1).+wsmb_dis);

			% 	Imw_yz=trapz(yz_dis(:,1),D_mat_mw(:,1).*yz_dis(:,1));
			% 	Imw_yz=Imw_yz.*2
			% ミス、対称だから慣性乗積はゼロ


%zz
zz_dis(:,1)=S_mat(:,1);
	%Imw_zz=trapz(S_mat(:,1),D_mat_mw(:,1).*(zz_dis(:,1)).^2);
Imw_zz=sum(zz_dis(:,1).^2.*mat_mass(:,1));		
Imw_zz=Imw_zz.*2

%xx
xx_dis(:,1)=sqrt(S_mat(:,1).^2.+(S_mat_swang(:,1).+wsmb_dis).^2);
Imw_xx=sum(xx_dis(:,1).^2.*mat_mass(:,1));
Imw_xx=Imw_xx.*2

			% %パイロット慣性モーメント計算
			% pw_dis=0.121;	%input("重心からパイロット重心までの距離(m)=");
			% p_mass=66;	%input("パイロットの体重(kg)=");

			% Ipw_yy=p_mass*pw_dis^2;
			% Ipw_xx=Ipw_yy

			% %ペラ慣性モーメントを計算
			% r_mass=1.2;
			% rw_dis=1.75;

			% Irw_yy=rw_dis^2*r_mass;
			% Irw_zz=Irw_yy


%各慣性モーメントを加算する
%ペラ
Ine_per=calc_ine(1.4,1.797,0,0.805)
%上駆動部
Ine_upd=calc_ine(0.65,1.399,0,0.799)
%下駆動部
Ine_lwd=calc_ine(1.16,0.919,0,-0.351)
%前タイヤ
Ine_ft=calc_ine(0.3,0.866,0,-0.757)
%後ろタイヤ
Ine_bt=calc_ine(0.3,-0.583,0,-0.757)
%フレーム
Ine_fl=calc_ine(3.95,-0.92,0,0.381)
%椅子とパイロットとフェアリング
Ine_pa=calc_ine(67.9,0.118,0,-0.387)
%テールパイプ
Ine_tl=calc_ine(1.4,-2.955,0,0.798)
%水平尾翼
Ine_ht=calc_ine(0.75,-4.513,0,0.909)
%垂直尾翼
Ine_vt=calc_ine(0.75,-5.079,0,0.849)

%全慣性モーメントに加算する

Ine(1,1)=Ine(1,1).+Imw_xx;
Ine(2,2)=Ine(2,2).+Imw_yy;
Ine(3,3)=Ine(3,3).+Imw_zz;

Ine=Ine+Ine_per+Ine_upd+Ine_lwd+Ine_ft+Ine_bt+Ine_fl+Ine_pa+Ine_tl+Ine_ht+Ine_vt