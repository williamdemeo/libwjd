function write_matrix(M, filename, mode)
if nargin < 3,
  mode = 'a';
end;
fid = fopen(filename,mode);
n=size(M,1);
for i=1:n,
  fprintf(fid,'%10.6f ',M(i,:));  fprintf(fid,'\n');
end;
fclose(fid);


