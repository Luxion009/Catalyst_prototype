%���^�̍��W�_�ǂݍ���
now_work=pwd;

read_d="FX76MP149_e66";
cd (read_d)

fpa=fopen('FX76MP149.dat');
fanamea=fgetl(fpa);#���^���ǂݍ���
foila=(fscanf(fpa,'%f',[2,Inf]))';#�Q�s������œǂݍ��񂾌�]�u
fpb=fopen('e66.dat');
fanameb=fgetl(fpb);#���^���ǂݍ���
foilb=(fscanf(fpb,'%f',[2,Inf]))';#�Q�s������œǂݍ��񂾌�]�u


front_p_a=0;
do
	front_p_a+=1;
until(foila(front_p_a,1)<=min(foila(:,1)))

front_p_b=0;
do
	front_p_b+=1;
until(foilb(front_p_b,1)<=min(foilb(:,1)))

cd (now_work)