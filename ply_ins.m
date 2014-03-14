%‚O“x‚ÌˆêüÏ‘w‚ğ‘}“ü
function ply_data=ply_ins(ply_data,insp,C_T,ply_end)

plyscount=rows(ply_data);
ply_len=ply_data(rows(ply_data),6);
ply_end=ply_end;

for n=2:6
	ply_buf=ply_data(insp:plyscount,n);
	ply_data(insp+1:plyscount+1,n)=ply_buf;

endfor

% cpe=rows(ply_data)+1;
% ply_data(cpe,2)=24;
% ply_data(cpe,3)=90;
% ply_data(cpe,4)=90;
% ply_data(cpe,5)=0;
% ply_data(cpe,6)=ply_len;


ply_data(insp,2)=40;
ply_data(insp,3)=90;
ply_data(insp,4)=0;
ply_data(insp,5)=0;
ply_data(insp,6)=ply_end;

for j=2:rows(ply_data)
	ply_data(j,1)=ply_data(j-1,1)+C_T(lookup(C_T(:,1),ply_data(j,2)),2)*2;
end


endfunction