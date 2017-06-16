from os.path import join
import json

from cogent import LoadTree
from numpy.testing import assert_array_less

import gapped
from data import get_data_dir
from . import get_aln

def test_joint_reconstruction():
    data_dir = get_data_dir()
    with open(join(data_dir, 'small_cnfgtr.json')) as lf_in:
        flat_cnfgtr = json.load(lf_in)
    model = lambda: gapped.CNFGTR(optimise_motif_probs=True, model_gaps=True)
    cnfgtr = gapped.inflate_likelihood_function(flat_cnfgtr, model)
    aln = get_aln(join(data_dir, 'small_aln.fasta'), 
                  filter_gaps=False, codon_position=-1)[:99]
    cnfgtr.setAlignment(aln)

    anc_aln = gapped.joint(cnfgtr)

    def prob_for_col(col):
        p = cnfgtr.getMotifProbs()[col['root']]
        p *= cnfgtr.getPsubForEdge('BterF3')[col['root']][col['BterF3']]
        p *= cnfgtr.getPsubForEdge('OsmaF4')[col['root']][col['OsmaF4']]
        p *= cnfgtr.getPsubForEdge('twoBeesF2')[col['root']][col['twoBeesF2']]
        p *= cnfgtr.getPsubForEdge('AmelF2')[col['twoBeesF2']][col['AmelF2']]
        p *= cnfgtr.getPsubForEdge('AdorF1')[col['twoBeesF2']][col['AdorF1']]
        return p

    flipper = {'A':'C', 'C':'G', 'G':'T', 'T':'A'}
    for i in range(0,99,3):
        col = anc_aln[i:i+3].todict()
        p = prob_for_col(col)
        if col['twoBeesF2'] == '---':
            col['twoBeesF2'] = 'CAC'
        else:
            col['twoBeesF2'] = \
                col['twoBeesF2'][:2] + flipper[col['twoBeesF2'][2]]
        pm = prob_for_col(col)
        assert_array_less(pm, p)
