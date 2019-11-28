#General

The earlier scripts (0-1) are mainly to read in the data from the JSON files.

Several participant numbers were replaced,
due to the mess discussed below. 
1.5-justify_replaced data reads that data and shows
why it was replaced

You can run from script 2-create_final_dfs.R to see most of the 
relevant details about subsequent data exclusions (where participant
number was not replaced)

Have added a column 'practice' - logical denoting whether
the trial is in a practice block or not. Practice block
includes the 100 AR practice trials they did and block 
1 of the stop task (184 trials).



2.5-check_balance checks how our counterbalance is going. 




#Testing messiness
One day we ran out of JSON files (at p number 50) and did not have the script to
generate new ones, so we needed some way to test the rest of our participnats that day.
Student ran participant 51 as a copy of p3 to match the counterbalance, but we soon 
realised this method not ideal because then 51 would see the trials in exact same randomized
order as p3. Thus for the rest of the day we replaced some participant numbers for participants who had
data issues - see Results/REPLACED

4 counterbalance orders based on key assignment/ condition order:
a,b,c,d

Runs in order a,b,c,d from files 1-51. However at 52 it starts
from 'b' and keeps going. This is because we had to re-run the script
to generate JSON files after maxing out at 50. The script to re-generate
JSON files after we ran out always started from 'a' regardless of participant number.
As mentioned above participant 51 used a re-used file which is why the balance swaps at 52