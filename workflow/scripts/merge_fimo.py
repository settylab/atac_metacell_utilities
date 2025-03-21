## adapted from https://gist.github.com/Mr-Milk/22c412454beb020141b041a125a2deaa
## Joins individually computed FIMO results to 
from pathlib import Path
import pandas as pd
import numpy as np

result_dir = Path('../../data/fimo_result/')
result_file = []
for f in result_dir.glob('individual_results/MA*'):
    result_file.append(f/'fimo.tsv')
    
cols = ['motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'strand',
       'score', 'p-value', 'q-value', 'matched_sequence']
data = pd.concat([pd.read_csv(r, sep="\t").iloc[0:-3, :] for r in result_file])[cols] #read in everything but the last 3 lines of each file and concatenate
data.to_csv(f'{result_dir.name}fimo.tsv', sep = '\t', index = False)