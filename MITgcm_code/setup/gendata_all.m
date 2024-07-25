clear all, close all

% First get parameters
params = gendata_params();
f = params.f;
deltaT = params.deltaT;
n_periods = params.n_periods;
k = params.k;
om = params.om;

flags.write_k_dependent = true;
flags.write_om_dependent = true;

for i = 1:length(om) % forcing frequency index
    om_prefix = sprintf('om%0.8f_',om(i)); % File prefix for om
    flags.write_om_dependent = true;

        % Set up run directories for every run
        rname = sprintf('run_om%0.8f',om(i));
        disp(rname)
        rdir = fullfile('..','runs',rname);
        if ~exist(rdir,'dir'); mkdir(rdir); end

        % Generate run data
        params = gendata(om(i), k(i), rdir, flags);

        % Template files have been set up with FIELD = PLACEHOLDER for fields that need
        % to be written. Make a cell array of substitutions with entries that
        % look like: {filename, {field1, string; ... fieldN, string}, prefix}
        other_subs = {};
        om_subs = {};
        k_subs = {};

        
        % Things we need to write once per forcing frequency
        if flags.write_om_dependent
            om_subs = {'../code/templates/SIZE.h',...
                     {'sNx', sprintf('%d', params.nxc/params.np(1));
                       'sNy', sprintf('%d', params.nyc/params.np(2));
                         'Nr', sprintf('%d', params.nzc)},...
                         om_prefix;
			'../input/templates/data.obcs',...
	                {'OB_Ieast', sprintf('%d*-1', params.nyc);
        	         'OB_Iwest', sprintf('%d*1', params.nyc);
                	 'OB_Jnorth', sprintf('%d*-1', params.nxc);
                	 'OB_Jsouth', sprintf('%d*1', params.nxc)},...
               		 om_prefix;
                    '../input/templates/data',...
                       {'deltaT', sprintf('%.1f',deltaT);
                        'nTimeSteps',  sprintf('%d',ceil(n_periods*2*pi/om(i)/deltaT));
                        'monitorFreq', sprintf('%.1f',2*pi/om(i)/8)},... 
                         om_prefix;
                        '../code/templates/rbcs_fields_load.F',...
                        {'om', sprintf('%0.8f',om(i));
                        'kf', sprintf('%0.8f',k(i));
                        },...
                        om_prefix;
			 '../input/templates/data.diagnostics',...
                       {'frequency(1)', sprintf('%.1f',-2*pi/om(i)/8);
                        'timePhase(1)', sprintf('%.1f',-2*pi/om(i)/8);
                        'frequency(2)', sprintf('%.1f',-2*pi/om(i)/8);
                        'timePhase(2)', sprintf('%.1f',-2*pi/om(i)/8);
                        'frequency(3)', sprintf('%.1f',-2*pi/om(i)/8);
                        'timePhase(3)', sprintf('%.1f',-2*pi/om(i)/8);
                        'frequency(4)', sprintf('%.1f',-2*pi/om(i)/8);
                        'timePhase(4)', sprintf('%.1f',-2*pi/om(i)/8)},...
                       om_prefix};

        end

        

        substitutions = cat(1,om_subs);

        for nf = 1:size(substitutions,1)
            ftxt = fileread(substitutions{nf,1}); % load template file text
            for ns = 1:size(substitutions{nf,2},1)
                expr = sprintf('%s *= *PLACEHOLDER',substitutions{nf,2}{ns,1}); % create regexp
                expr = strrep(expr,'(','\('); % Replace parens for regexp function
                expr = strrep(expr,')','\)');
                [i1,i2] = regexp(ftxt,expr,'start','end'); % find placeholder

                % update placeholder text with substitution text
                ftxt = cat(2, ...
                           ftxt(1:i1-1), ...
                           strrep(ftxt(i1:i2),'PLACEHOLDER',substitutions{nf,2}{ns,2}),...
                           ftxt(i2+1:end));
            end
            % Write new text to file named with prefix
            [fdir,fname,fext] = fileparts(substitutions{nf,1});
            dir_out = strrep(fdir,'templates','generated');
            newfile = fullfile(dir_out,sprintf('%s%s%s',substitutions{nf,3},fname,fext));
            fid = fopen(newfile,'w');
            fprintf(fid,'%s\n',ftxt);
            fclose(fid);
        end % substitutions in template files

        % run shell script to link files
        cmd=sprintf('kwics_setup.sh %s %s %s',om_prefix,'--build');
        system(cmd);

        flags.write_om_dependent = false; % We've written data for this om
    
end % loop over om
