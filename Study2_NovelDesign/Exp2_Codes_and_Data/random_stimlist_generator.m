    function [list1, list2] = random_stimlist_generator(n,number_of_trials)

    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % ================= created by Tom Salomon, January 2015 ==================
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    % This function is used to create a set of random binary comparisons to be
    % later used in binary_ranking function.
    
    % ---------------- Input ----------------
    % n - number of stimuli.
    % number_of_trials - desired number of unique comparisons.
    
    % ---------------- Output ---------------
    % [list1,list2] - two number_of_trials sized vectors, in which for every i
    % [list1(i),list2(i)] is a unique pair; list1(i)<list2(i).
    
    % These comparisons are all unique - the same comparison may not return
    % more than once. In addition the function will try to make sure each
    % stimulus appears an even amount of time.
    % In case it is not possible to present each stimulus an even amount of
    % times, (e.g. in case the user asked 61 trials with only 60
    % stimuli), the function will not fail, but it will print a warning.
    
    % ================================
    % IMPORTANT NOTE I:
    % Every stimulus will be presented 2*number_of_trials/n times.
    % Therefore, for equal presentation of each stimulus, number_of_trials
    % must be a multiple of n/2.
    %  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % e.g if n=60, number_of_trials can be 30,60,90,120....
    % ================================
    
    % for testing:
    % n=60;
    % number_of_trials=300;
    % % -----------------------------------------------------------------------
    
    % essential for randomization
    rng('shuffle');
    
    time_start=GetSecs;
    
    % This part will print a warning if the function cannot generate lists with
    % equal number of presentations for each stimulus.
    if rem(2*number_of_trials,n)~=0
        fprintf(2,'<strong>random_stimlist_generator WARNING:</strong> not all stimuli will be presented an even number of times!\n');
        fprintf(2,'The function will not fail, however - this may cause problems further down the road.\n \n');
        fprintf(2,'For even presentation, make sure number_of_trials variable is a multiple of n/2\n');
        fprintf(2,'e.g if n=60, number_of_trials should be 30,60,90,120...\n');
    end
    
    % Create empty variables to be filled.
    presentations_count=zeros(n,2); % n by 2 matrix in which first column indicate the stimulus, and the 2nd column the number of appearances.
    presentations_count(:,1)=(1:n);
    list1=zeros(1,number_of_trials);
    list2=zeros(1,number_of_trials);
    trial=1;
    restart=0;
    
    while trial<=number_of_trials
        %check how much time passed since code started
        time_passed=GetSecs-time_start;
        
        % This part will print a warning if the function cannot generate lists
        % due to too many comparisons and too few stimuli.
        if number_of_trials>0.5*n*(n-1)
            fprintf(2,'\n<strong>random_stimlist_generator ERROR:</strong> pairing is impossible!\n');
            fprintf(2,'number of trials is too large in comparison to number of stimuli \n');
            fprintf(2,'maximun number of comparisons = n(n-1)/2 \n');
            list1=[];
            list2=[];
            break;
        end
        
        % locate all stimuli which were presented the least amount of times.
        % From these stimuli, 2 will be randomly chosen and paired
        stimuli=Shuffle(presentations_count(presentations_count(:,2)==min(presentations_count(:,2)),1));
        
        % In case only one stimulus is left to choose from, select this stimulus
        % and pair it with another stimulus which was presented one time more
        % than this stimulus. Normally, this part is in only active if 'n' is an
        % odd number.
        if length(stimuli)<2
            add_stimulus=stimuli;
            stimuli=Shuffle(presentations_count(presentations_count(:,2)==(min(presentations_count(:,2))+1),1));
            stimuli(end+1)=stimuli(1);
            stimuli(1)=add_stimulus;
        end
        
        % randomly select two stimuli (by now stimuli variable is shuffled)
        stimuli=sort(stimuli(1:2));
        
        % make sure that this pair was not selected already
        if sum(list2(list1==stimuli(1))==stimuli(2))==0
            presentations_count(presentations_count(:,1)==stimuli(1),2)=presentations_count(presentations_count(:,1)==stimuli(1),2)+1;
            presentations_count(presentations_count(:,1)==stimuli(2),2)=presentations_count(presentations_count(:,1)==stimuli(2),2)+1;
            
            list1(trial)=min(stimuli);
            list2(trial)=max(stimuli);
            trial=trial+1; % pairing is OK, continue to next trial
            
            % Untie condition: In case only one possible pair is left, and this pair is illegal
            % (already appeared)"untie" this case. go back one selection,
            % remove the previous pairing and re-pair
        else
            if length(presentations_count(presentations_count(:,2)==min(presentations_count(:,2)),1))<=2
                trial=trial-1;
                presentations_count(presentations_count(:,1)==list1(trial),2)=presentations_count(presentations_count(:,1)==list1(trial),2)-1;
                presentations_count(presentations_count(:,1)==list2(trial),2)=presentations_count(presentations_count(:,1)==list2(trial),2)-1;
                list1(trial)=0;
                list2(trial)=0;
                restart=restart+1; %
                
            end
            % Restart condition: In the rare cases the "untie" procedure did
            % not work (the untie loop repeated 10 times), restart the entire code.
            % For very exuastive procedures (very small n/number_of_trials ratio)
            % consider increasing the restart condition bar (the "if restart==X" below).
            if restart==10 || time_passed>1
                presentations_count=zeros(n,2);
                presentations_count(:,1)=(1:n);
                list1=zeros(1,number_of_trials);
                list2=zeros(1,number_of_trials);
                trial=1;
                restart=0;
                time_start=GetSecs;
                
            end
        end
    end % end of trial loop.
end


