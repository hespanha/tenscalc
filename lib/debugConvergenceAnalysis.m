function debugConvergenceAnalysis(debugInfo,status,iter,time,z,nu,lambda,dZ,dNu,dLambda,G,F)
% debugConvergence(debugInfo,status,iter,time,z,nu,lambda,dZ,dNu,dLambda,G,F)
%
% Takes as inputs all the outputs from an ipmPD solver and provides
% suggestions to improve the convergence.
%
%
% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.

% study scaling of z
    
    zBlocks=[cellfun(@(x)prod(size(x)),debugInfo.P1optimizationVariables),...
             cellfun(@(x)prod(size(x)),debugInfo.P2optimizationVariables)];
    k1=cellfun(@(x)strcmp(x.type,'iszero'),debugInfo.P1constraints);
    k2=cellfun(@(x)strcmp(x.type,'iszero'),debugInfo.P2constraints);
    GBlocks=[cellfun(@(x)prod(size(x.objs{1})),debugInfo.P1constraints(k1)),...
             cellfun(@(x)prod(size(x.objs{1})),debugInfo.P2constraints(k2))];
    k1=cellfun(@(x)strcmp(x.type,'ispositive'),debugInfo.P1constraints);
    k2=cellfun(@(x)strcmp(x.type,'ispositive'),debugInfo.P2constraints);
    FBlocks=[cellfun(@(x)prod(size(x.objs{1})),debugInfo.P1constraints(k1)),...
             cellfun(@(x)prod(size(x.objs{1})),debugInfo.P2constraints(k2))];
             
    fig=55;

    figure(fig);set(gcf,'name','Dispersion');fig=fig+1;
    checkScaling(z(:,1),zBlocks,'z = optim. vars: initial',...
                 z(:,end),zBlocks,'z: final',...
                 nu(:,1),GBlocks,'nu - eq. constr. multipliers:initial',...
                 nu(:,end),GBlocks,'nu: final',...
                 lambda(:,1),FBlocks,'lambda - ineq constr. mult.: initial',...
                 lambda(:,end),FBlocks,'lambda: final');
        
    fprintf('Problem: alpha -> 0\n');
    fprintf('Possibly solution: large dispersion in values of optimization varibles (or in \n');
    fprintf('                   dual variables), causing very large steps in some variables\n');

    fprintf('Rules regarding dispersion of optimization variables (z):\n');
    fprintf('. max ration for initial & final values should not be much larger than 1000.\n');
    fprintf('  If this is not the case, rescale optimization variables accordingly.\n');
    fprintf('  (replacing the variable by X times itself in the objectives and constraints,\n');
    fprintf('   divides the variable by X).\n')

    fprintf('Rules regarding dispersion of equality constraints multipliers (nu):\n');
    fprintf('. largest final values should not be much larger than 1000.\n');
    fprintf('  If this is not the case, rescale the equality constraints appropriately\n');
    fprintf('  (multiplying both side of equality by X, divides nu by X).\n');
    fprintf('. smallest final value close to zero means that there are dependent equalities.\n');
    
    fprintf('Rules regarding dispersion of inequality constraints multipliers (lambda):\n');
    fprintf('. largest final values should not be much larger than 1000.\n');
    fprintf('  If this is not the case, rescale the equality constraints appropriately\n');
    fprintf('  (multiplying both side of inequalty by X>0, divides lambda by X).\n');

    fprintf('Problem: gradient and/equality do not -> 0\n');
    fprintf('Possibly solution: ?\n');

end 

function checkScaling(varargin)
    
    clf
    nPlots=length(varargin)/3;
    s=1;
    for i=1:3:length(varargin)
        data=varargin{i};
        blocks=varargin{i+1};
        name=varargin{i+2};
        if length(data)~=sum(blocks)
            size(data)
            blocks
            error('size of data does not match block lengths\n');
        end
        blockID=nan(length(data));
        blockNdx=nan(length(data));
        k=1;
        for i=1:length(blocks)
            blockID(k:k+blocks(i)-1)=i;
            blockNdx(k:k+blocks(i)-1)=1:blocks(i);
            k=k+blocks(i);
        end
        
        subplot(2,nPlots,s);s=s+nPlots;
        semilogy(abs(data));
        grid on
        title(sprintf('%s - magnitudes',name))
        subplot(2,nPlots,s);s=s-nPlots+1;
        semilogy(sort(abs(data)),'.')
        grid on
        title(sprintf('%s - sorted magnitudes',name))
        
        fprintf(sprintf('Dispersion of %s\n',name));

        % find extreme entries
        knz=find(abs(data)>0);
        [m,km]=min(abs(data(knz)));
        [M,kM]=max(abs(data(knz)));
        kSmall=knz(km);
        kLarge=knz(kM);
        fprintf('   . max ratio = %10.2g -> between entries %3d (entry %3d in block %2d = %10.2g) and %3d (%3d entry in block %2d = %10.2g)\n',M/m,...
                kSmall,blockNdx(kSmall),blockID(kSmall),data(kSmall),...
                kLarge,blockNdx(kLarge),blockID(kLarge),data(kLarge));
        
        if 0
            % find gap between entries
            [ds,ks]=sort(abs(data));
            r=ds(2:end)./ds(1:end-1);
            r(r==inf)=0;
            [maxGapFactor,km]=max(r);
            kSmall=ks(km);
            kLarge=ks(km+1);
            
            fprintf('   . max gap   = %10.2g -> between entries %3d (entry %3d in block %2d = %10.2g) and %3d (%3d entry in block %2d = %10.2g)\n',maxGapFactor,...
                    kSmall,blockNdx(kSmall),blockID(kSmall),data(kSmall),...
                    kLarge,blockNdx(kLarge),blockID(kLarge),data(kLarge));
        end
    end
end
    