function isSensitivity=variableIndices(u,optimizationVariables,whereVariables,sensitivityVariables)
% find indices of sensitivity variables in u

    isSensitivity=false(length(u),1);

    for i=1:length(sensitivityVariables)
        if strcmp(type(sensitivityVariables{i}),'subsref')
            S=parameters(sensitivityVariables{i});
            sensitivityVariables{i}=Tcalculus(operands(sensitivityVariables{i}));
            ind=1:prod(size(sensitivityVariables{i}));
            ind=reshape(ind,msize(sensitivityVariables{i}));
            ind=subsref(ind,S);
        else
            ind=-1;
        end
        if ~strcmp(type(sensitivityVariables{i}),'variable')
            sensitivityVariables{i}
            error('all sensitivityVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(sensitivityVariables{i}));
        end
        found=false;
        for j=1:length(optimizationVariables)
            if isequal(sensitivityVariables{i},optimizationVariables{j})
                if ind<0
                    fprintf('   sensitivityVariable %s: %d values\n',sensitivityVariables{i}.name,length(whereVariables{j}));
                    isSensitivity(whereVariables{j})=true;
                else
                    fprintf('   sensitivityVariable %s[...]: %d values\n',sensitivityVariables{i}.name,length(whereVariables{j}(ind)));
                    isSensitivity(whereVariables{j}(ind))=true;
                end
                found=true;
                break;
            end
        end
        if ~found
            error('sensitivityVariable %s is not an optimizationVariable\n',name(sensitivityVariables{i}));
        end
    end
    if any(isSensitivity)
        fprintf('      %d sensitivity variables\n',sum(isSensitivity));
    end
end