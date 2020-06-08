import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
input_trimmed = sys.argv[1]
input_untrimmed = sys.argv[2]
output = sys.argv[3]

def get_df(input):
    all_reads = []
    with open(input) as f:
      for line in f:
        if line.count('\t') == 0: #context line or '====' or empty line
          if line.count('context') == 1: #context line
            context = line.strip().split(' context ')[0]
            read = line.strip().split(' context ')[1][1:3]
        else:
          if line.strip().split('\t')[0] != 'position': #exclude header line
            to_add = [0.0 if x=='' else float(x) for x in line.strip().split('\t')]
            to_add += [context, read]
            all_reads.append(to_add)
    df = pd.DataFrame(all_reads, columns = ['position', 'count_methylated', 'count_unmethylated', 'perc_methylation', 'coverage', 'context', 'read'])
    return(df)

df = get_df(input_untrimmed)
fig, axs = plt.subplots(3,2, sharex=True, sharey=False, figsize=(8,11))
fig.suptitle('Methylation bias by position\nuntrimmed                      trimmed')

axs[0,0].plot( df[(df['context']=='CpG') & (df['read']=='R1')]['position'], df[(df['context']=='CpG') & (df['read']=='R1')]['perc_methylation'], label='Read 1' )
axs[0,0].plot( df[(df['context']=='CpG') & (df['read']=='R2')]['position'], df[(df['context']=='CpG') & (df['read']=='R2')]['perc_methylation'], label='Read 2' )
axs[0,0].set_title('CpG')
axs[0,0].set_ylabel('% methylation')
axs[0,0].set_yticks(np.arange(40, 100, 5))
axs[0,0].legend(loc="upper right")

axs[1,0].plot( df[(df['context']=='CHG') & (df['read']=='R1')]['position'], df[(df['context']=='CHG') & (df['read']=='R1')]['perc_methylation'], label='Read 1' )
axs[1,0].plot( df[(df['context']=='CHG') & (df['read']=='R2')]['position'], df[(df['context']=='CHG') & (df['read']=='R2')]['perc_methylation'], label='Read 2' )
axs[1,0].set_title('CHG')
axs[1,0].set_ylabel('% methylation')
axs[1,0].set_yticks(np.arange(0, 2, 0.5))
axs[1,0].legend(loc="upper right")

axs[2,0].plot( df[(df['context']=='CHH') & (df['read']=='R1')]['position'], df[(df['context']=='CHH') & (df['read']=='R1')]['perc_methylation'], label='Read 1' )
axs[2,0].plot( df[(df['context']=='CHH') & (df['read']=='R2')]['position'], df[(df['context']=='CHH') & (df['read']=='R2')]['perc_methylation'], label='Read 2' )
axs[2,0].set_title('CHH')
axs[2,0].set_xlabel('Position in read (bp)')
axs[2,0].set_ylabel('% methylation')
axs[2,0].set_yticks(np.arange(0, 2, 0.5))
axs[2,0].legend(loc="upper right")

df = get_df(input_trimmed)
axs[0,1].plot( df[(df['context']=='CpG') & (df['read']=='R1')]['position'], df[(df['context']=='CpG') & (df['read']=='R1')]['perc_methylation'], label='Read 1' )
axs[0,1].plot( df[(df['context']=='CpG') & (df['read']=='R2')]['position'], df[(df['context']=='CpG') & (df['read']=='R2')]['perc_methylation'], label='Read 2' )
axs[0,1].set_title('CpG')
axs[0,1].set_ylabel('% methylation')
axs[0,1].set_yticks(np.arange(40, 100, 5))
axs[0,1].legend(loc="upper right")

axs[1,1].plot( df[(df['context']=='CHG') & (df['read']=='R1')]['position'], df[(df['context']=='CHG') & (df['read']=='R1')]['perc_methylation'], label='Read 1' )
axs[1,1].plot( df[(df['context']=='CHG') & (df['read']=='R2')]['position'], df[(df['context']=='CHG') & (df['read']=='R2')]['perc_methylation'], label='Read 2' )
axs[1,1].set_title('CHG')
axs[1,1].set_ylabel('% methylation')
axs[1,1].set_yticks(np.arange(0, 2, 0.5))
axs[1,1].legend(loc="upper right")

axs[2,1].plot( df[(df['context']=='CHH') & (df['read']=='R1')]['position'], df[(df['context']=='CHH') & (df['read']=='R1')]['perc_methylation'], label='Read 1' )
axs[2,1].plot( df[(df['context']=='CHH') & (df['read']=='R2')]['position'], df[(df['context']=='CHH') & (df['read']=='R2')]['perc_methylation'], label='Read 2' )
axs[2,1].set_title('CHH')
axs[2,1].set_xlabel('Position in read (bp)')
axs[2,1].set_ylabel('% methylation')
axs[2,1].set_yticks(np.arange(0, 2, 0.5))
axs[2,1].legend(loc="upper right")

plt.xticks(np.arange(0,72, 5))
plt.savefig(output)
