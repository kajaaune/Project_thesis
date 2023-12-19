clear

dir_cm = 'D:/Project_thesis/Whole_data_matlab2'; 
plotx2east;
plotzIntoPlane;

lattice_parameters_post = [8.08 8.08 8.08]; % Post cycle 8.17 before
lattice_parameters = [8.17 8.17 8.17];
angles = [90 90 90] * degree;

cs = crystalSymmetry('m-3m', lattice_parameters, angles, 'mineral',...
    'LNMO', 'X||a*', 'Z||c*'); % Note the X||a* alignment

cs_post = crystalSymmetry('m-3m', lattice_parameters, angles, 'mineral',...
    'LNMO', 'X||a*', 'Z||c*'); % Note the X||a* alignment

% Load data
% Pristine samples
% ebsd1 = EBSD.load(fullfile(dir_cm,...
%         "xmap1_whole_v2.ang"), cs,'setting 2','convertEuler2SpatialReferenceFrame');
% ebsd2 = EBSD.load(fullfile(dir_cm,...
%         "xmap2_whole_v2.ang"), cs,'setting 2','convertEuler2SpatialReferenceFrame');
% ebsd3 = EBSD.load(fullfile(dir_cm,...
%         "xmap3_whole_v2.ang"), cs,'setting 2','convertEuler2SpatialReferenceFrame');

% Post_cycle samples
ebsd1p = EBSD.load(fullfile(dir_cm,...
        "xmap1_post_whole.ang"), cs_post,'setting 2','convertEuler2SpatialReferenceFrame');
ebsd2p = EBSD.load(fullfile(dir_cm,...
        "xmap2_post_whole.ang"), cs_post,'setting 2','convertEuler2SpatialReferenceFrame');
ebsd3p = EBSD.load(fullfile(dir_cm,...
        "xmap_post3_v2.ang"), cs_post,'setting 2','convertEuler2SpatialReferenceFrame');
% ebsd1.scanUnit = 'nm';
% ebsd2.scanUnit = 'nm';
% ebsd3.scanUnit = 'nm';
ebsd1p.scanUnit = 'nm';
ebsd2p.scanUnit = 'nm';
ebsd3p.scanUnit = 'nm';

%%
% Index shifts: (use when choosing grains)
% Tilt 1 is standard
% Tilt 2 has x-shift -1 and y-shift -1
% Tilt 3 has x-shift -4 and y-shift -7

%%
% Calculate grains and big grains
% [grains1,big_g1] = make_grains(ebsd1);
% [grains2,big_g2] = make_grains(ebsd2);
% [grains3,big_g3] = make_grains(ebsd3);

[grains1p,big_g1p] = make_grains(ebsd1p);
[grains2p,big_g2p] = make_grains(ebsd2p);
[grains3p,big_g3p] = make_grains(ebsd3p);
%%
% Plot the big grains
% plot_grains(big_g1,ebsd1);
% plot_grains(big_g2,ebsd2);
% plot_grains(big_g3,ebsd3);

plot_grains(big_g1p,ebsd1p);
plot_grains(big_g2p,ebsd2p);
plot_grains(big_g3p,ebsd3p);

%%
% table = table_grain_properties(big_g1);
tablep = table_grain_properties(big_g1p);

%%
% First check the elements of tilt 2 to tilt 1
% table = remove_grains(table,big_g2,big_g3);
tablep = remove_grains(tablep,big_g2p,big_g3p);
%%
% The grain id of all the grains we want
% grain_id1 = [];
% grain_id2 = [];
% grain_id3 = [];
% for index = 1:length(table.x)
%     grain_id1 = [grain_id1 big_g1(table.x(index),table.y(index)).id];
%     grain_id2 = [grain_id2 big_g2(table.x(index)-1,table.y(index)-1).id];
%     grain_id3 = [grain_id3 big_g3(table.x(index)-4,table.y(index)-7).id];
% end

grain_id1p = [];
grain_id2p = [];
grain_id3p = [];
for index = 1:length(tablep.x)
    grain_id1p = [grain_id1p big_g1p(tablep.x(index),tablep.y(index)).id];
    grain_id2p = [grain_id2p big_g2p(tablep.x(index)-1,tablep.y(index)-1).id];
    grain_id3p = [grain_id3p big_g3p(tablep.x(index)-4,tablep.y(index)-7).id];
end

%%
% Make tables for tilt 2 and 3
% table2 = table_grain_properties(grains2(grain_id2));
% table3 = table_grain_properties(grains3(grain_id3));

table2p = table_grain_properties(grains2p(grain_id2p));
table3p = table_grain_properties(grains3p(grain_id3p));
%%
% Plot the common grains for each tilt
% figure
% plot(grains1(grain_id1))
% figure
% plot(grains2(grain_id2))
% figure
% plot(grains3(grain_id3))

%%
% Write to txt file for processing in python
%%{
% writetable(table,'D:/Matlab_scripts/Matlab_ori_data2/tilt1_oris_v2.txt','Delimiter',' ')  
% 
% writetable(table2,'D:/Matlab_scripts/Matlab_ori_data2/tilt2_oris_v2.txt','Delimiter',' ')  
% 
% writetable(table3,'D:/Matlab_scripts/Matlab_ori_data2/tilt3_oris_v2.txt','Delimiter',' ')  


writetable(tablep,'D:/Matlab_scripts/Matlab_ori_data2/tilt1p_oris_v3.txt','Delimiter',' ')  

writetable(table2p,'D:/Matlab_scripts/Matlab_ori_data2/tilt2p_oris_v3.txt','Delimiter',' ')  

writetable(table3p,'D:/Matlab_scripts/Matlab_ori_data2/tilt3p_oris_v3.txt','Delimiter',' ')  
%%}

%%
disp(height(tablep));
disp(height(table2p));
disp(height(table3p));

%%

% Function for removing grains not present in all three tilts
function table = remove_grains(table,big2,big3)
remove_indices = [];
for index = 1:length(table.x)
   tilt2_grain = big2(table.x(index)-1,table.y(index)-1);
   tilt3_grain = big3(table.x(index)-4,table.y(index)-7);
   if isempty(tilt2_grain.id) == true
      remove_indices = [remove_indices index];
   else

   end
   if isempty(tilt3_grain.id) == true
      remove_indices = [remove_indices index];
   end
end
remove_indices = unique(remove_indices);
table(remove_indices,:) = [];
disp(length(remove_indices));
end

function T = table_grain_properties(big_grains)
    xy = big_grains.centroid;
    x = xy(:,1);
    y = xy(:,2);
    e1 = big_grains.meanOrientation.phi1;
    e2 = big_grains.meanOrientation.Phi;
    e3 = big_grains.meanOrientation.phi2;
    T = table(x,y,e1,e2,e3);
end

function [grains,big_grains] = make_grains(data)
    [grains, data.grainId, data.mis2mean] = calcGrains(data, 'angle',...
        3*degree); % 3 looks good
    grains = smooth(grains,50);
    grains_LNMO = grains('LNMO');
    big_grains = grains_LNMO(grains_LNMO.grainSize > 80); % exclude smallest grains % 80 seems fine
end

function plotted = plot_grains(big_grains,ebsd)
    directions = {xvector, yvector, zvector,vector3d(1,1,1)};
    ipfkey = ipfTSLKey(ebsd.CS); 

    
    ipfkey.inversePoleFigureDirection = directions{3};
    omcolor = ipfkey.orientation2color(ebsd.orientations);
    %figure
    plot(ipfkey)
    figure
    plot(ebsd, omcolor)
    hold on
    legend('off')
    %plot(big_grains,cScube*0.5);
    plot(big_grains.boundary,'lineWidth',2);
    %text(big_grains,int2str(big_grains.id));
    hold off
    plotted = 0;
end
