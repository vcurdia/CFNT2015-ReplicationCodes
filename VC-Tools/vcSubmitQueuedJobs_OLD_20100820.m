function JobOut = vcSubmitQueuedJobs(nJobs,JobFcnHandle,nJobArgOut,JobOptions,...
    ListPathDependencies,nMaxWorkers,ShowOutput,LogFileName,OutputPreamble)

% vcSubmitQueuedJobs
%
% Submit parallel jobs that might exceed the maximum number of workers
% available.
%
% Usage:
%   vcSubmitQueuedJobs(nJobs,JobFcnHandle,nJobArgOut,JobOptions,...
%       ListPathDependencies)
%   vcSubmitQueuedJobs(...,nMaxWorkers)
%   vcSubmitQueuedJobs...,nMaxWorkers,ShowOutput)
%   vcSubmitQueuedJobs(...,nMaxWorkers,ShowOutput,LogFileName)
%   vcSubmitQueuedJobs(...,nMaxWorkers,ShowOutput,LogFileName,OutputPreamble)
%
% Mandatory Inputs:
%   
%   nJobs
%   Number of Jobs to run
%
%   JobFcnHandle
%   Handle to the function called by each worker.
%
%   nJobArgOut
%   Number of output arguments in each job.
%   
%   JobOptions
%   Cell array with options required by the function called by the worker.
%   It needs to contain one element per job, unless no inputs are
%   necessary, in which case, {} needs to be used.
%
%   ListPathDependencies
%   Cell array with the list of path dependencies needed in function called
%   by workers.
%
%   nMaxWorkers
%   [Optional] Maximum number of simultaneous workers. Default: 4.
%
%   ShowOutput
%   If set to 1 all output of each worker is shown. If LogFileName is not
%   specified, then it is shown in the caller command history. If
%   LogFileName is specified then it is stored in log.
%   Default: 1
%
%   LogFileName
%   If specified, then all output is saved in a ascii file.
%
%   OutputPreamble
%   If specified this should be a string to introduce in a fprintf command
%   before showing the output of each worker. It must receive as input the
%   worker number.
%
% .........................................................................
%
% Created: February 11, 2009 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2009-2011 by Vasco Curdia

%% ------------------------------------------------------------------------

if ~exist('nMaxWorkers','var'), nMaxWorkers = 4; end
if ~exist('ShowOutput','var'), ShowOutput = 1; end
if ~exist('LogFileName','var')||isempty(LogFileName),Save2Log=0;else Save2Log=1;end
if isempty(JobOptions),for j=1:nJobs,JobOptions{j} = {};end,end
nCompleted = 0;
nSubmitted = 0;
j = 0;
jobsRunning = false(nJobs,1);
JobOut = cell(nJobs,1);
while nCompleted<nJobs
    while (nSubmitted<nMaxWorkers) && (j<nJobs)
        % submit new jobs
        j = j+1;
        nSubmitted = nSubmitted+1;
        job{j} = createSgeJob(0);
        set(job{j},'PathDependencies',ListPathDependencies)
        TaskID{j} = createTask(job{j},JobFcnHandle,nJobArgOut,JobOptions{j});
        if ShowOutput
            set(TaskID{j}, 'CaptureCommandWindowOutput', true);
        end
        submit(job{j})
        jobsRunning(j) = true;
    end
    pause(0.01)
    listJobsRunning = find(jobsRunning);
    for jj=1:length(listJobsRunning)
        if strcmp(get(job{listJobsRunning(jj)},'State'),'finished')
            jobsRunning(listJobsRunning(jj)) = false;
            nSubmitted = nSubmitted-1;
            nCompleted = nCompleted+1;
        end
    end
end
for j=1:nJobs
    if ShowOutput
        if Save2Log
            jobLogName = sprintf('%s%.0f.log',LogFileName,j);
            fid = fopen(jobLogName,'wt');
            if exist('OutputPreamble','var')
                fprintf(fid,OutputPreamble,j);
            end
            fprintf(fid,strrep(get(TaskID{j},'CommandWindowOutput'),'%','%%'));
            fclose(fid);
        else
            if exist('OutputPreamble','var')
                fprintf(OutputPreamble,j);
            end
            fprintf(strrep(get(TaskID{j},'CommandWindowOutput'),'%','%%'));
        end
    end
    if ~isempty(get(TaskID{j},'errormessage'))
        fprintf(strrep(get(TaskID{j},'erroridentifier'),'%','%%'));fprintf('\n');
        fprintf(strrep(get(TaskID{j},'errormessage'),'%','%%'));fprintf('\n');
    end
    JobOut{j} = getAllOutputArguments(job{j});
    destroy(job{j})
end
JobOut = cat(1,JobOut{:});
