import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import seaborn as sns
sns.set_style('ticks')
sns.set_context('paper')

from glob import glob
import sys

from RiboGraphViz import RGV
sys.path.append('/Users/hwayment/das/github/DegScore')
from DegScore import DegScore

keyword = sys.argv[1]
input_files = sorted(glob(keyword + '*RUNNING_BEST.txt'))
if len( input_files ) == 0:
    print( 'Could not find any input files beginning with keyword',keyword)
    exit(0)

data = pd.DataFrame()
for fil in input_files:
    cond = fil.split('.')[0]
    identifier = fil.split('.')[-3].split('-')[-1]
    tmp = pd.read_csv(fil,delimiter='\t')
    tmp['condition'] = cond
    tmp['id'] = identifier
    data = data.append(tmp,ignore_index=True)

condition_list = list(data.condition.unique())
palette = sns.color_palette('tab10', len(condition_list))

# Plot statistics of run

plt.figure(figsize=(10,4))
nrows, ncols = 2,3

plt.subplot(nrows,ncols,1)
for job in data.id.unique():
    tmp = data.loc[data.id==job]
    condition = tmp.condition.iloc[0]
    plt.plot(tmp.DegScore.values, color=palette[condition_list.index(condition)])

plt.ylabel('DegScore')
plt.xlabel('Checkpoint')

plt.subplot(nrows,ncols,2)
for job in data.id.unique():
    tmp = data.loc[data.id==job]
    condition = tmp.condition.iloc[0]
    plt.plot(tmp.AUP.values, color=palette[condition_list.index(condition)])

plt.ylabel('AUP')
plt.xlabel('Checkpoint')

plt.subplot(nrows, ncols, 3)
sns.scatterplot(x='DegScore',y='dG(MFE)', data=data, hue='condition',palette=palette)
plt.legend(bbox_to_anchor=(1,1), frameon=False)

plt.subplot(nrows, ncols, 4)
sns.scatterplot(x='DegScore',y='CAI', data=data, hue='condition',palette=palette)
plt.legend([], frameon=False)

plt.subplot(nrows, ncols, 5)
sns.scatterplot(x='DegScore',y='AUP_init14', data=data, hue='condition',palette=palette)
plt.legend([], frameon=False)
plt.ylim([0,1])

plt.subplot(nrows, ncols, 6)
sns.scatterplot(x='DegScore',y='MLD', data=data, hue='condition',palette=palette)
plt.legend([], frameon=False)

plt.tight_layout()

plt.savefig('%s_output.pdf' % keyword, bbox_inches='tight')
print('Saved %s_output.pdf' % keyword)

n_rows = len(condition_list)

plt.figure(figsize=(10,15))

for cond_ind, condition in enumerate(condition_list):

    tmp = data.loc[data.condition==condition]
    tmp = tmp.sort_values('DegScore')
    tmp = tmp.drop_duplicates('id')

    for i in range(3):
        plt.subplot(n_rows, 3, cond_ind*3+i+1)

        print('Generating RiboGraphViz for %s, %d/3' % (condition, i+1))

        struct = tmp.iloc[i]['MFE Structure']
        seq = tmp.iloc[i]['full_sequence']

        mdl = DegScore(seq, structure=struct)
        degscore_vector = mdl.degscore_by_position

        rgv_object = RGV(struct)
        rgv_object.draw(c=degscore_vector)

        plt.title("%s\nRun ID: %s\nDegscore: %d" % (condition, tmp.iloc[i]['id'], tmp.iloc[i]['DegScore']))

plt.savefig('%s_best_3_from_each_condition.pdf' % keyword, bbox_inches='tight')
print('Saved %s_best_3_from_each_condition.pdf' % keyword)
