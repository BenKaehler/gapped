from numpy import prod, argmax
import cogent

def joint(lf, locus=None):
    aln = lf.getParamValue('alignment', locus=locus)
    collen = len(list(lf.model.getAlphabet())[0])
    anc_aln = None
    for i in range(0,len(aln),collen):
        col = reconstruct(lf, aln[i:i+collen], lf.tree, locus=locus)
        anc_aln = col if anc_aln is None else anc_aln + col
    return anc_aln

def reconstruct(lf, aln, tree, locus=None):
    if tree.isTip():
        y = tree.Name
        P = lf.getPsubForEdge(y, locus=locus)
        j = str(aln.NamedSeqs[y])
        if 'N' in j:
            js = lf.model.getAlphabet().resolveAmbiguity(j)
            L = lambda i: max(P[i,j] for j in js)
            def C(i):
                j = argmax([P[i,j] for j in js])
                return js[j]
        else:
            L = lambda i: P[i,j]
            C = lambda i: j
        tree.C = C
        return L
    
    Ls = [reconstruct(lf, aln, c, locus=locus) for c in tree.Children]
    alphabet = list(lf.model.getAlphabet())
    calcedLs = {j:prod([L(j) for L in Ls]) for j in alphabet}
    
    if tree.isRoot():
        pi = lf.getMotifProbs()
        j = argmax([pi[j]*calcedLs[j] for j in alphabet])
        tree.anc = alphabet[j]
        result = [(tree.Name, tree.anc)]
        for child in tree.Children:
            _get_anc(child, result)
        return cogent.LoadSeqs(data=result)
        
    P = lf.getPsubForEdge(tree.Name, locus=locus)
    L = lambda i: max(P[i,j]*calcedLs[j] for j in alphabet)
    def C(i):
        j = argmax([P[i,j]*calcedLs[j] for j in alphabet])
        return alphabet[j]
    tree.C = C
    return L

def _get_anc(tree, result):
    tree.anc = tree.C(tree.Parent.anc)
    result.append((tree.Name, tree.anc))
    for child in tree.Children:
        _get_anc(child, result)
