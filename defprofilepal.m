%���^�쐬��⏕����

%cl*c���`
for i=1:(wld_c+we_div)
	c_cl(i,1)=2*ed_Hex_gamma(i,1)/ed_Hex_comp_vel(i,1);
endfor


do

	%�嗃���E���C�m���Y����`
	lim_mre=input("�嗃���E���C�m���Y��=");

	%���E���C�m���Y������ŊO���������`
	chord(we_div,1)=lim_mre*nu/d_vel;
	disp("�ŊO���[����");
	disp(chord(we_div,1));

	%���������t���p
	ang_cw=input("�������L���}�p��");
	mt_ang=ang_cw+ed_Hex_i_ang(1,1)

	%�������悭������`
	chord(1:dw_divn,1)=input("������������(m)=");

	%�Ȃ�������
	chord(dw_divn*2,1)=input("����������(m)=");

	%�O���[
	chord(dw_divn*3,1)=input("�O���[������(m)=");

	%���ŊO���[
	chord(dw_divn*4,1)=input("���ŊO���[����(m)=");
	
	disp("�������X��");
	for i=2:w_div_c
		i
		chord_alpha(i,1)=(chord(dw_divn*i,1)-chord(dw_divn*(i-1),1))/divw_span(i,1);
		disp(chord_alpha(i,1));

	endfor

	%���������z���쐬
	disp("generating chord");
	for i=2:w_div_c
		chord_d(i,1)=(chord(dw_divn*(i),1)-chord(dw_divn*(i-1),1))/dw_divn;
		chord((i-1)*dw_divn+1,1)=chord((i-1)*dw_divn,1)+chord_d(i,1);
		for j=2:dw_divn
			chord((i-1)*dw_divn+j,1)=chord((i-1)*dw_divn+j-1,1)+chord_d(i,1);
		endfor
	endfor

	%���C�m���Y�����z�쐬
	disp("generating re");
	chord_re(:,1)=chord(:,1).*d_vel/nu;

	%Cl���z�쐬
	disp("generating cl");
	chord_cl(1:we_div,1)=c_cl(1:we_div,1)./chord(1:we_div,1);


	figure(10);
	plot(nang_line(1:we_div),chord(1:we_div),"-");

	for i=1:w_div_c
		printf("%d��\n",i)
		disp("Re��");
		disp(chord_re(dw_divn*i,1));
		disp("cl");
		disp(chord_cl(dw_divn*i,1));
		disp("�U���}���p")
		disp(ed_Hex_i_ang(dw_divn*i,1));
		disp("�L���}�p")
		disp(mt_ang-ed_Hex_i_ang(dw_divn*i,1))
		disp("---------------------")
	endfor

	chord_sig=input("1:ok 0:return ---");

until(chord_sig==1)

%�e���[�K�v�f�[�^���Z�o