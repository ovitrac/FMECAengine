% http://depth-first.com/articles/2007/01/24/thirty-two-free-chemistry-databases/
% http://www.chemspider.com/inchi.asmx
% http://www.chemspider.com/inchi.asmx/InChIToCSID?inchi=#{inchi}

% http://www.chemspider.com/AboutServices.aspx
% Chemspider Services (with/without security token)
createClassFromWsdl('http://www.chemspider.com/inchi.asmx?WSDL')  % http://www.chemspider.com/inchi.asmx
createClassFromWsdl('http://www.chemspider.com/Search.asmx?WSDL') % http://www.chemspider.com/Search.asmx
createClassFromWsdl('http://www.chemspider.com/MassSpecAPI.asmx?WSDL')  % http://www.chemspider.com/MassSpecAPI.asmx

%% Usage
token = '30c079d0-cedb-42a7-be40-6e8f9f2c0d75';
obj = Search;
csid=SimpleSearch(obj,'BHT',token);
nfo = GetCompoundInfo(obj,csid.int,token);
% img
base64string=GetCompoundThumbnail(obj,csid.int,token);
encoder = org.apache.commons.codec.binary.Base64;
img = encoder.decode(uint8(base64string));
fid = fopen(fullfile(tempdir,'tst.png'),'w'); fwrite(fid,img,'int8'); fclose(fid);

% advanced example)
csid = SimpleSearch(obj,'irganox 1076',token);
nfo = GetCompoundInfo(obj,csid.int{1},token);


