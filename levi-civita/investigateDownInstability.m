%A script which investigates the sudden instability of the down orbits in
%the interpolated regime in the vicinity of the maximum and minimum. Run
%fixedRadiusAnalysisInterp first.

figure;
hold on;

%indexmin = mini-1;
%indexmax = mini+1;
numPts = 20;
%thetasinvest = linspace(endubdthetas(indexmin),endubdthetas(indexmax),numPts);
thetasinvest = linspace(0.804,0.806,numPts);

%We interpolate the down points integrated backwards in order to find down
%points which more exactly correspond with up points.
enddbdintinvest = spline(enddbdthetas,enddbdsorted,thetasinvest);
enddbdintpolarinvest = cart2m2polar(enddbdintinvest, mu);
enddbdintthetasinvest = enddbdintpolarinvest(2,:);

intdownlclistinvest = zeros(4,numel(theta));
enddfdint = zeros(4,numel(theta));
solintdnfdEndinvest = {};

c = startCaching(c);

%We call getFromCache for the sole purpose of caching the equations of
%motion.
[~,c] = getFromCache(c,'ca.eqnsHandle',@defaultHandle,c);

%We now integrate the interpolated down intercepts back through the origin
%and to the detection circle once more.
parfor i = 1:numel(enddbdintthetasinvest)
    
    disp(i)
    
    %Now, we must convert back to Levi-Civita form:
    [~,intdownlc] = standard2lc(0,enddbdintinvest(:,i),mu);
    
    solintdnfdinvest=integ(linspace(0,2*endtime,1000),intdownlc,c);
    solintdnfdEndinvest{i} = solintdnfdinvest;

end

c = stopCaching(c);

disp('Integration complete.')

for i = 1:numel(enddbdintthetasinvest)
    timesdnfd = linspace(solintdnfdEndinvest{i}.x(1),solintdnfdEndinvest{i}.x(end),10000);
    [~,trajdnfdpts] = lc2standard(timesdnfd,deval(timesdnfd,solintdnfdEndinvest{i}),...
                                        mu);
    cplot(trajdnfdpts,c)
end

