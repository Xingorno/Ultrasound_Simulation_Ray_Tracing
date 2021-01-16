fileID = 'prelog_rf.txt';
prelog_rf_temp = textread(fileID,'%f');
prelog_rf = reshape(prelog_rf_temp, [2048, 465])';
figure(1)
imshow(prelog_rf)
%%
fileID = 'postlog_rf.txt';
postlog_rf_temp = textread(fileID,'%f');
postlog_rf = reshape(postlog_rf_temp, [2048, 465])';
figure(2)
imshow(postlog_rf)

%%

fileID = 'Simulated_US.txt';
US_temp = textread(fileID,'%f');
US = reshape(US_temp, [500, 400])';
figure(8)
imshow(US)

