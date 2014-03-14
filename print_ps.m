%}–Êo—Í

%—ƒŒ^Ä“Ç‚İ‚İ
fp=fopen('use_foila.txt','r');
use_foil_ps{1}=fgetl(fp);
fclose(fp);

fp=fopen('use_foilb.txt','r');
use_foil_ps{2}=fgetl(fp);
fclose(fp);

fp=fopen('use_foilc.txt','r');
use_foil_ps{3}=fgetl(fp);
fclose(fp);

fp=fopen('use_foild.txt','r');
use_foil_ps{4}=fgetl(fp);
fclose(fp);

disp("input profiledata_ps");

for i=1:4
	%-----—ƒŒ^‚ğ“Ç
	i
	fp=fopen(strcat(use_foil_ps{i},'.dat'));
	buff=fgetl(fp);
	foil_ps{i}=(fscanf(fp,'%f',[2,350]))';
	foil_ps{i}=(foil_ps{i});
	fclose(fp);
	fs(i)=0;

	%-----—ƒ‚Ì‘O‰‚ğ’Tõ
	do	
		fs(i)=fs(i)+1;
	until(foil_ps{i}(fs(i),1)<=min(foil_ps{i}(:,1)))

end
	
%-----—ƒŒ^ƒf[ƒ^‚ğ®—
%-----x•ûŒü•ªŠ„”
n_xdiv_l=80;
n_xdiv_t=60;

disp("generating xy_foil_ps");

x_foil_ps(1:n_xdiv_t+n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([n_xdiv_t+n_xdiv_l-1:-1:1]').^2;	
x_foil_ps(n_xdiv_t+n_xdiv_l,:)=0;
x_foil_ps(n_xdiv_t+n_xdiv_l+1:2*n_xdiv_t+2*n_xdiv_l-1,1)=1./(n_xdiv_t+n_xdiv_l-1)^2.*([1:n_xdiv_t+n_xdiv_l-1]').^2;

for i=1:4
	y_foil_ps(1:n_xdiv_l+n_xdiv_t,i)=interp1(foil_ps{i}(1:fs(i),1),foil_ps{i}(1:fs(i),2),x_foil_ps(1:n_xdiv_l+n_xdiv_t,1),"extrap");
	y_foil_ps(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,i)=interp1(foil_ps{i}(fs(i):rows(foil_ps{i}),1),foil_ps{i}(fs(i):rows(foil_ps{i}),2),x_foil_ps(n_xdiv_l+n_xdiv_t:2*n_xdiv_t+2*n_xdiv_l-1,1),"extrap");
	y_foil_ps(n_xdiv_l+n_xdiv_t,i)=0;
end

%o—Í
disp("output A wing")

%A—ƒ

end_pp=print_ex(dwing_chord(1,1),dwing_aoa(1,1),dwing_spos(1,1),dwing_sdia(1,:),wing1_cutr,[x_foil_ps(:,1) y_foil_ps(:,1)],'C',1,1,pj_name,now_work);

%B—ƒ
disp("B wing")
for i=1:divnum(2,1)+1
	i
	buf_yfoil=y_foil_ps(:,1).*print_mix(2,i).+y_foil_ps(:,2).*(1-print_mix(2,i));
	buf_cut=wing1_cutr.*print_mix(2,i).+wing2_cutr.*(1-print_mix(2,i));
	end_pp=print_ex(print_chord(2,i),print_aoa(2,i),print_spos(2,i),dwing_sdia(2,:),buf_cut,[x_foil_ps(:,1) buf_yfoil(:,1)],'N',i,print_mix(2,i),pj_name,now_work);
	clear buf_yfoil buf_cut
endfor

%C—ƒ
disp("C wing")
for i=1:divnum(3,1)+1
	i
	buf_yfoil=y_foil_ps(:,2).*print_mix(3,i).+y_foil_ps(:,3).*(1-print_mix(3,i));
	buf_cut=wing2_cutr.*print_mix(3,i).+wing3_cutr.*(1-print_mix(3,i));
	end_pp=print_ex(print_chord(3,i),print_aoa(3,i),print_spos(3,i),dwing_sdia(3,:),buf_cut,[x_foil_ps(:,1) buf_yfoil(:,1)],'G',i,print_mix(3,i),pj_name,now_work);
	clear buf_yfoil buf_cut
endfor

%D—ƒ
disp("D wing")
for i=1:divnum(4,1)+1
	i
	buf_yfoil=y_foil_ps(:,3).*print_mix(4,i).+y_foil_ps(:,4).*(1-print_mix(4,i));
	buf_cut=wing3_cutr.*print_mix(4,i).+wing4_cutr.*(1-print_mix(4,i));
	end_pp=print_ex(print_chord(4,i),print_aoa(4,i),print_spos(4,i),dwing_sdia(4,:),buf_cut,[x_foil_ps(:,1) buf_yfoil(:,1)],'FSG',i,print_mix(4,i),pj_name,now_work);
	clear buf_yfoil buf_cut
endfor

disp("E wing")
for i=1:divnum(5,1)+1
	i
	buf_yfoil=y_foil_ps(:,4);
	buf_cut=0
	end_pp=print_ex(print_chord(5,i),print_aoa(5,i),print_spos(5,i),dwing_sdia(5,:),buf_cut,[x_foil_ps(:,1) buf_yfoil(:,1)],'SSG',i,0,pj_name,now_work);
	clear buf_yfoil buf_cut
endfor

