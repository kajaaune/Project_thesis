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

directions = {xvector, yvector, zvector};

% Pristine samples
ebsd1 = EBSD.load(fullfile(dir_cm,...
        "xmap1_whole_v2.ang"), cs,'setting 2','convertEuler2SpatialReferenceFrame');
ebsd2 = EBSD.load(fullfile(dir_cm,...
        "xmap2_whole_v2.ang"), cs,'setting 2','convertEuler2SpatialReferenceFrame');
ebsd3 = EBSD.load(fullfile(dir_cm,...
        "xmap3_whole_v2.ang"), cs,'setting 2','convertEuler2SpatialReferenceFrame');
% Post_cycle samples
ebsd1p = EBSD.load(fullfile(dir_cm,...
        "xmap1_post_whole.ang"), cs_post,'setting 2','convertEuler2SpatialReferenceFrame');
ebsd2p = EBSD.load(fullfile(dir_cm,...
        "xmap2_post_whole.ang"), cs_post,'setting 2','convertEuler2SpatialReferenceFrame');
ebsd3p = EBSD.load(fullfile(dir_cm,...
        "xmap_post3_v2.ang"), cs_post,'setting 2','convertEuler2SpatialReferenceFrame');
ebsd1.scanUnit = 'nm';
ebsd2.scanUnit = 'nm';
ebsd3.scanUnit = 'nm';

ebsd1p.scanUnit = 'nm';
ebsd2p.scanUnit = 'nm';
ebsd3p.scanUnit = 'nm';
%%
% Make big grains and identify twinning
[grains1,big_g1,gb_1,twins1,twinning1] = make_grains(ebsd1);
[grains2,big_g2,gb_2,twins2,twinning2] = make_grains(ebsd2);
[grains3,big_g3,gb_3,twins3,twinning3] = make_grains(ebsd3);

[grains1p,big_g1p,gb_1p,twins1p,twinning1p] = make_grains(ebsd1p);
[grains2p,big_g2p,gb_2p,twins2p,twinning2p] = make_grains(ebsd2p);
[grains3p,big_g3p,gb_3p,twins3p,twinning3p] = make_grains(ebsd3p);

%%
% Plotting the twin boundaries
plot_twins(ebsd1,grains1,big_g1,twins1)
plot_twins(ebsd2,grains2,big_g2,twins2)
plot_twins(ebsd3,grains3,big_g3,twins3)

plot_twins(ebsd1p,grains1p,big_g1p,twins1p)
plot_twins(ebsd2p,grains2p,big_g2p,twins2p)
plot_twins(ebsd3p,grains3p,big_g3p,twins3p)

%%
% Find rotation axis
% Twinning help https://github.com/mtex-toolbox/mtex/discussions/975
gb = gb_1p;
ind = gb.misorientation.angle>59*degree & gb.misorientation.angle<61*degree;
mori = gb.misorientation(ind);
mori_mean = mean(mori,'robust');
r2m = round2Miller(mori_mean);
%%
% prisine
% twinning = orientation.map(Miller({2 -1 2},cs),Miller({0 0 1},cs),Miller({1 0 -1},cs),Miller({-1 -1 0},cs));
% cycled
twinning = orientation.map(Miller({1 2 1},cs),Miller({-1 -1 2},cs),Miller({1 0 -1},cs),Miller({-1 1 0},cs));

% Roation axis
rot_ax = round(twinning.axis);

% Rotational angle
twin_angle = twinning.angle / degree;
%%
% Line profile for long grains
path = ['C:/Users/Kaja Eggen Aune/OneDrive - NTNU/Bilder uni/Project_thesis_images/' ...
    'Results/Misorientations_within_tilt/Twin_lines/post_line_twins_1.png'];
ebsd_l = plot_line(ebsd1p,grains1p,path);
%%
% Function for calculating grains, finding big grains, and twinning
function [grains,big_grains,gB_LNMO,twinBoundary,twinning] = make_grains(data)
    [grains, data.grainId, data.mis2mean] = calcGrains(data, 'angle',...
        3*degree); % 3 looks good
    grains = smooth(grains,50);
    
    grains_LNMO = grains('LNMO');
    big_grains = grains_LNMO(grains_LNMO.grainSize > 80); % exclude smallest grains % 80 seems fine
    CS = grains.CS;
    gB = big_grains.boundary;
    gB_LNMO = gB('LNMO','LNMO');
    % Twinning
    % Pristine tilt 1
    twinning = orientation.map(Miller(1,1,-1,CS),Miller(-1,-1,-1,CS),...
    Miller(1,-1,0,CS,'uvw'),Miller(-1,0,1,CS,'uvw'));  
    % Cycled tilt 1
    % twinning = orientation.map(Miller(1,-1,1,CS),Miller(1,1,1,CS),...
    % Miller(1,0,-1,CS,'uvw'),Miller(-1,-1,0,CS,'uvw'));  
    isTwinning = angle(gB_LNMO.misorientation,twinning) < 4*degree;
    twinBoundary = gB_LNMO(isTwinning);
    

end

% Function for plotting the twin boundaries
function plt = plot_twins(ebsd,grains,big_grains,twins)
    directions = {xvector, yvector, zvector};
    ipfkey = ipfTSLKey(ebsd.CS);
    ipfkey.inversePoleFigureDirection = directions{3};
    omcolor = ipfkey.orientation2color(ebsd.orientations);
    % plot the twinning boundaries
    figure
    plot(grains,grains.meanOrientation,omcolor)
    %plot(ebsd('indexed'),ebsd('indexed').orientations)
    hold on
    %plot(gB_MgMg,angle(gB_LNMO.misorientation,twinning),'linewidth',4)
    plot(twins,'linecolor','w','linewidth',4)%,'displayName','twin boundary')
    hold off
    plt = 0;
end

% Function for calculating and plotting the long grains, line and plot
% misorientation angle and orientation gradient.
function ebsd_l = plot_line(ebsd,grains,path)
    directions = {xvector, yvector, zvector};
    ipfkey = ipfTSLKey(ebsd.CS);
    ipfkey.inversePoleFigureDirection = directions{3};
    omcolor = ipfkey.orientation2color(ebsd.orientations);
    % Plot the grains
    % Long grain 1
    % grains_l = [4442,4083,4019];
    % lineSec =  [2687 2939; 2318 4714];
    
    % Long grain 2
    % grains_l = [2669,2595,2518,1908,1981];%,2027,1545,1231]; 
    % lineSec =  [1557 4205; 823, 2778];
    
    % Long grain 3
    % grains_l = [5928,5930,5649,4972,4671];
    % lineSec =  [3612 2587; 3066 3958];

    % Cycled
    grains_l = [2028,1634,1436];
    lineSec =  [1672 1033; 1308 1054];
    
    ebsd_line = spatialProfile(ebsd,lineSec);
    vec = vector3d(lineSec(3)-lineSec(1),lineSec(4)-lineSec(2),0);   
    test = inv(ebsd_line.orientations)*vec;

    % draw the line into the plot
    figure
    plot(grains(grains_l),grains(grains_l).meanOrientation,omcolor)  
    line(lineSec(:,1),lineSec(:,2),'linewidth',2)
    
    figure
    plot(ebsd_line.x,...
    angle(ebsd_line(1).orientations,ebsd_line.orientations)/degree,'linewidth',3)

    % misorientation gradient
    hold all
    plot(0.5*(ebsd_line.x(1:end-1)+ebsd_line.x(2:end)),...
      angle(ebsd_line(1:end-1).orientations,ebsd_line(2:end).orientations)/degree,'linewidth',3)
    hold off
    ax = gca; 
    ax.FontSize = 24;
    ylim(ax,[0 65])
    xlabel('x','FontSize',26); ylabel('misorientation angle in degree','FontSize',26)

    legend('to reference orientation','orientation gradient','FontSize',26,'Location','northeast')

    ebsd_l = test;%ebsd_line;
end