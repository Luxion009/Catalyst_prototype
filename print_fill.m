%�}�ʏo�͗p�Ƀf�[�^��ۊǂ���
%�����͂��ׂ�mm�\�L

disp("data fill")
%�������[������
for i=1:5
	dwing_chord(i,1)=real_chord(i,1)*1000;
endfor

%�}�p
for i=1:5
	dwing_aoa(i,1)=f_mtang(i*50,1);
endfor

%���ʒu
for i=1:5
	dwing_spos(i,1)=spar_pos(i*50,1);
endfor

%���a
%(:,1)�ő�a(:,2)�ŏ��a

dwing_sdia=[111.4 110.7;103.0 101.5;82.6 81.3;61.9 61.1;39.7 39.7];

%������
wing1_cut=[-222.62 -7.99 73.85;-124.35 -0.3 67.29;119.44 -7.46 77.92;206.11 -17.57 66.6;277.35 -28.83 47.64];
wing2_cut=[-204.59 -8.57 64.31;-115.13 -0.81 67.27;103.27 -6.08 67.91;187 -15.63 68.63;261.46 -27.66 50.22];
wing3_cut=[-162.87 -2.34 53.49;-89.73 0.6 56.9;84.3 -4.15 56.25;153.37 -10.12 54.41;214.63 -17.35 40.96];
wing4_cut=[-119.7 -0.63 39.21;-66.36 0.98 43.52;67.24 -3.54 45.44;123.25 -8.19 42.97;172.82 -13.68 32.79];

disp("rotate cutpoint")
%�������ʒu���}���p�Ȃ��ɕϊ�
cut_t=3.69052/180*pi;
rot_sq=[cos(cut_t) -sin(cut_t);sin(cut_t) cos(cut_t)];
wing1_cutr=wing1_cut;
for i=1:5
	wing1_cutr(i,1:2)=(rot_sq*wing1_cutr(i,1:2)')';
endfor

wing2_cutr=wing2_cut;
for i=1:5
	wing2_cutr(i,1:2)=(rot_sq*(wing2_cutr(i,1:2))')';
endfor

wing3_cutr=wing3_cut;
for i=1:5
	wing3_cutr(i,1:2)=(rot_sq*(wing3_cutr(i,1:2))')';
endfor

wing4_cutr=wing4_cut;
for i=1:5
	wing4_cutr(i,1:2)=(rot_sq*(wing4_cutr(i,1:2))')';
endfor


disp("filling dwing data")
%B���f�[�^�ۊ�
%���u�Ԑ�
divnum(2,1)=14;

for i=1:divnum(2,1)+1
	print_chord(2,i)=dwing_chord(1,1)-((dwing_chord(1,1)-dwing_chord(2,1))/divnum(2,1))*(i-1);
endfor

for i=1:divnum(2,1)+1
	print_aoa(2,i)=dwing_aoa(2,1);
endfor

for i=1:divnum(2,1)+1
	print_spos(2,i)=dwing_spos(1,1)-((dwing_spos(1,1)-dwing_spos(2,1))/divnum(2,1))*(i-1);
endfor

%������
%��������1
for i=1:divnum(2,1)+1
	print_mix(2,i)=1-(1/divnum(2,1))*(i-1);
endfor


%C���f�[�^�ۊ�
divnum(3,1)=14;

for i=1:divnum(3,1)+1
	print_chord(3,i)=dwing_chord(2,1)-((dwing_chord(2,1)-dwing_chord(3,1))/divnum(3,1))*(i-1);
endfor

for i=1:divnum(3,1)+1
	print_aoa(3,i)=dwing_aoa(3,1);
endfor

for i=1:divnum(3,1)+1
	print_spos(3,i)=dwing_spos(2,1)-((dwing_spos(2,1)-dwing_spos(3,1))/divnum(3,1))*(i-1);
endfor

for i=1:divnum(3,1)+1
	print_mix(3,i)=1-(1/divnum(3,1))*(i-1);
endfor

%D���f�[�^�ۊ�
divnum(4,1)=13;

for i=1:divnum(4,1)+1
	print_chord(4,i)=dwing_chord(3,1)-((dwing_chord(3,1)-dwing_chord(4,1))/divnum(4,1))*(i-1);
endfor

for i=1:divnum(4,1)+1
	print_aoa(4,i)=dwing_aoa(4,1);
endfor

for i=1:divnum(4,1)+1
	print_spos(4,i)=dwing_spos(3,1)-((dwing_spos(3,1)-dwing_spos(4,1))/divnum(4,1))*(i-1);
endfor

for i=1:divnum(4,1)+1
	print_mix(4,i)=1-(1/divnum(4,1))*(i-1);
endfor

%E���f�[�^�ۊ�
divnum(5,1)=9;

for i=1:divnum(5,1)+1
	print_chord(5,i)=dwing_chord(4,1)-((dwing_chord(4,1)-dwing_chord(5,1))/divnum(5,1))*(i-1);
endfor

for i=1:divnum(5,1)+1
	print_aoa(5,i)=dwing_aoa(4,1)-((dwing_aoa(4,1)-dwing_aoa(5,1))/divnum(5,1))*(i-1);
endfor

for i=1:divnum(5,1)+1
	print_spos(5,i)=dwing_spos(4,1)-((dwing_spos(4,1)-dwing_spos(5,1))/divnum(5,1))*(i-1);
endfor

for i=1:divnum(5,1)+1
	print_mix(5,i)=1-(1/divnum(5,1))*(i-1);
endfor

%������
ofc=strcat('ribdata_C','.txt');

fpc=fopen(ofc,'wt');
fprintf(fpc,'������\n')

fprintf(fpc,'���u�ԍ�:%d �R�[�h�F%f ���t���p:%f ���a:%f-%f ���ʒu:%f\n',1,dwing_chord(1,1),dwing_aoa(1,1),dwing_sdia(1,1),dwing_sdia(1,2),dwing_spos(1,1));


fclose(fpc);

%����
ofc=strcat('ribdata_N','.txt');

fpc=fopen(ofc,'wt');
fprintf(fpc,'����\n')
for i=1:divnum(2,1)+1
	fprintf(fpc,'���u�ԍ�:%d �R�[�h�F%f ���t���p:%f ���a:%f-%f ���ʒu:%f\n',i,print_chord(2,i),print_aoa(2,i),dwing_sdia(2,1),dwing_sdia(2,2),print_spos(2,i));
endfor

fclose(fpc);

%�O��
ofc=strcat('ribdata_G','.txt');

fpc=fopen(ofc,'wt');
fprintf(fpc,'�O��\n')
for i=1:divnum(3,1)+1
	fprintf(fpc,'���u�ԍ�:%d �R�[�h�F%f ���t���p:%f ���a:%f-%f ���ʒu:%f\n',i,print_chord(3,i),print_aoa(3,i),dwing_sdia(3,1),dwing_sdia(3,2),print_spos(3,i));
endfor

fclose(fpc);

%1�ŊO��
ofc=strcat('ribdata_FSG','.txt');

fpc=fopen(ofc,'wt');
fprintf(fpc,'���ŊO��\n')
for i=1:divnum(4,1)+1
	fprintf(fpc,'���u�ԍ�:%d �R�[�h�F%f ���t���p:%f ���a:%f-%f ���ʒu:%f\n',i,print_chord(4,i),print_aoa(4,i),dwing_sdia(4,1),dwing_sdia(4,2),print_spos(4,i));
endfor

fclose(fpc);

%2�ŊO��
ofc=strcat('ribdata_SSG','.txt');

fpc=fopen(ofc,'wt');
fprintf(fpc,'���ŊO��\n')
for i=1:divnum(5,1)+1
	fprintf(fpc,'���u�ԍ�:%d �R�[�h�F%f ���t���p:%f ���a:%f-%f ���ʒu:%f\n',i,print_chord(5,i),print_aoa(5,i),dwing_sdia(5,1),dwing_sdia(5,2),print_spos(5,i));
endfor

fclose(fpc);