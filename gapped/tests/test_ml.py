import json
import os
from nose.tools import assert_equal

from numpy.testing import assert_almost_equal
from cogent import LoadTree

from data import get_data_dir
from gapped import _ml
from gapped.nest import inflate_likelihood_function,\
    deflate_likelihood_function
import gapped
from . import get_aln

def test_GNC():
    with open(os.path.join(get_data_dir(), 'GNC.json')) as infile:
        flat_lf = json.load(infile)

    lf = inflate_likelihood_function(flat_lf, _ml.GNC)
    aln = get_aln(os.path.join(get_data_dir(),
        'ENSG00000100393.fasta.gz'), codon_position=-1)
    lf.setAlignment(aln)

    flat_again = deflate_likelihood_function(lf)

    assert_almost_equal(flat_lf['EN'].values(), flat_again['EN'].values(), 9)

def test_gapped_CNFGTR():
    aln = get_aln(os.path.join(get_data_dir(),
                  'ENSG00000100393.fasta.gz'), 
                  codon_position=-1, filter_gaps=False)
    tree = LoadTree(treestring='(Human,Mouse,Opossum);')
    doc = {'aln' : str(aln), 'tree' : str(tree)}
    cnfgtr_result = gapped.ml(doc, model='CNFGTR', model_gaps=True, 
                              omega_indep=False, indel_indep=False)
    model = lambda: gapped.CNFGTR(optimise_motif_probs=True, model_gaps=True)
    cnfgtr = gapped.inflate_likelihood_function(cnfgtr_result['lf'], model)

    pi = cnfgtr.getMotifProbsByNode()['root'].asarray()
    P = cnfgtr.getPsubForEdge('Human')
    assert_almost_equal(pi.dot(P), pi)

    omega = cnfgtr.getParamValue('omega')
    pi = cnfgtr.getMotifProbs()
    Q = cnfgtr.getRateMatrixForEdge('Human')
    cond_p = pi['CCG'] / sum(pi['CC'+c] for c in 'ACGT')
    ref_cell = Q['CCT']['CCG']/cond_p
    cond_p = pi['CCC'] / sum(pi['CC'+c] for c in 'ACGT')
    assert_almost_equal(Q['CCA']['CCC']/cond_p/ref_cell, 
                    cnfgtr.getParamValue('A/C'))
    assert_almost_equal(Q['---']['CCC']/pi['CCC']/ref_cell,
                    cnfgtr.getParamValue('indel'))
    R = Q.asarray()/pi.asarray()
    assert_almost_equal(R.T, R)
