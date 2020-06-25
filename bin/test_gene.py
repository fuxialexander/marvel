import numpy as np
from scipy.sparse import load_npz
import argparse

# argument
parser = argparse.ArgumentParser()
parser.add_argument("distance_breaks", default="",type=str,
                    help="Type of region, e.g. enhancer or promoter")
parser.add_argument("result_path", default="",type=str,
                    help="path of chunk testing results")
parser.add_argument("output_path", default="./", type=str, help="output directory")
args = parser.parse_args()


distance_breaks = np.array([
    0, 50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000,
    500000, 550000, 600000, 650000, 700000, 750000, 800000, 850000, 900000,
    950000, 1000000
])

distance_weight = np.array([
    0.132257992, 0.133288330, 0.135888400, 0.126786539, 0.107038924,
    0.085863692, 0.065999800, 0.051969109, 0.038726515, 0.029763539,
    0.022819252, 0.017063568, 0.013074640, 0.010261396, 0.007829280,
    0.006324146, 0.004990197, 0.004069675, 0.003359097, 0.002625910
])




class Results(object):
    """The results class. Containing labels, pvalues, degree of freedoms, motif
    coefficients and raw data matrices Xs.
    """

    def __init__(self,
                 filename=None,
                 override=False,
                 labels=None,
                 pvalues=None,
                 dofs=None,
                 coefs=None,
                 xs=None):
        """Results constructor.

        Arguments:
        - `filename`: Name of .npz file to load
        - `override`: Whether to overrive specific argument.
        - `labels`: region /group labels
        - `pvalues`: pvalues, len(c_arrays) * len(region/group)
        - `dofs`: degree of freedoms
        - `coefs`: motif coefficients
        - `xs`: raw data matrices
        """
        self._selection_pvalue = None
        self._qvalues = None
        self._selection_fdr = None

        if (filename is None and override is False):
            print("Need filename or data to load.")
            return -1

        if override:
            if labels:
                self._labels = labels
            if pvalues:
                self._pvalues = pvalues
            if dofs:
                self._dofs = dofs
            if coefs:
                self._coefs = coefs
            if xs:
                self._xs = xs

        else:
            with np.load(filename) as data:
                keys = data.keys()
                print(keys)
                key_rs = [i for i in keys if re.match('.*rs.*', i) is not None]
                key_ps = [i for i in keys if re.match('.*ps.*', i) is not None]
                key_ds = [i for i in keys if re.match('.*ds.*', i) is not None]
                key_cs = [i for i in keys if re.match('.*cs.*', i) is not None]
                key_xs = [i for i in keys if re.match('.*xs.*', i) is not None]
                if len(key_rs) > 0:
                    self._labels = data[key_rs[0]]
                if len(key_ps) > 0:
                    self._pvalues = data[key_ps[0]]
                if len(key_ds) > 0:
                    self._dofs = data[key_ds[0]]
                if len(key_cs) > 0:
                    self._coefs = data[key_cs[0]]
                if len(key_xs) > 0:
                    self._xs = data[key_xs[0]]

    def select_fdr(self, alpha=0.1):
        """Use BH FDR threshold to select results, output q-values and selected
        indices.

        Arguments:
        - `alpha`: FDR cut-off

        """
        self._selection_fdr, self._qvalues = fdrcorrection0(
            self._pvalues.min(1), alpha)
        self._selection_fdr = np.where(self._selection_fdr)[0]

    def select_pvalue(self, alpha=0.00005):
        """Use P-value threshold to select results, output selected indices.

        Arguments:
        - `alpha`: P-value cut-off

        """
        self._selection_pvalue = np.unique(np.where(self._pvalues < alpha)[0])

    def label_to_hg38(self, label_hg19_hg38_mapping_file, delim='\t'):
        """Produce labels in hg38. label_hg19_hg38_mapping_file should be a
        tab-delimited file with 1st column been hg19 coordinates and 2nd column
        been hg38 coordinates.

        Mapping file produced by:
        awk '{print $4,$1":"$2"-"$3}' ../hg38_encc_enhancer_null_4column.bed > label_hg19_hg38_mapping_enhancer.txt
        grep chr ../hg38_encc_enhancer_null_4column.bed.unmap| cut -f4 | awk '{print $0" Unmapped"}' >> label_hg19_hg38_mapping_enhancer.txt
        """
        mapping = {}
        with open(label_hg19_hg38_mapping_file, 'r') as f:
            for line in f:
                (key, val) = line.strip('\n').split(delim)
                mapping[key] = val

        self._labels_hg38 = np.array([
            mapping[i] if mapping.has_key(i) else "Unmapped"
            for i in self._labels
        ])

def get_pe(pe_pair_file):
    ep_dis = dict()
    ep_pair = dict()
    with open(pe_pair_file, 'r') as pairs:
        for pair in pairs:
            p, e, d = pair.strip("\n").split("\t")
            ep_dis[p + e] = int(d)
            if ep_pair.has_key(p):
                ep_pair[p].append(e)
            else:
                ep_pair[p] = []
                ep_pair[p].append(e)
    return ep_dis, ep_pair

def get_pg(pg_pair_file):
    pg_pair = dict()
    glist = []
    with open(pg_pair_file, 'r') as f:
        for line in f:
            loc, gid= line.strip('\n').split('\t')
            glist.append(gid)
            if pg_pair.has_key(gid):
                pg_pair[gid].append(loc)
            else:
                pg_pair[gid] = []
                pg_pair[gid].append(loc)
    return pg_pair, np.unique(np.array(glist))

enhancer_results = Results("enhancer_results.npz")

promoter_results = Results("promoter_results.npz")

enhancer_profiles = Results("enhancer_profiles.npz")

promoter_profiles = Results("promoter_profiles.npz")

ep_dis, ep_pair = get_pe("promoter_enhancer_pair.txt")

pg_pair, glist = get_pg("promoter_gene_pair.txt")

def get_weight(distance, breaks, weight):
    diff = breaks - abs(distance)
    if sum(diff > 0) > 0:
        return weight[np.where(diff > 0)[0][diff[diff > 0].argmin()] - 1]
    else:
        return 0

def get_distance(p, e, ep_dis):
    if ep_dis.has_key(p + e):
        return ep_dis[p + e]
    else:
        return 1000000 - 1

def get_gene(pg_pair, glist, ep_pair, sample_count, motif_count):
    xs = []
    for i, g in enumerate(glist):
        x = np.zeros((sample_count, motif_count))
        promoters = np.array(pg_pair[g])
        enhancers = [
            np.array(ep_pair[p]) for p in promoters if ep_pair.has_key(p)
        ]
        if len(enhancers) > 0:
            enhancers = np.unique(np.concatenate(enhancers))
        else:
            continue
        pids = np.array([
            np.where(promoter_results._labels == promoter)[0][0]
            for promoter in promoters
        ])
        for enhancer in enhancers:
            eid = np.where(enhancer_results._labels == enhancer)[0][0]
            weight = get_weight(
                get_distance(promoters[0], enhancer, ep_dis), distance_breaks,
                distance_weight)
            x += weight * enhancer_profiles[eid * sample_count:
                                      (eid + 1) * sample_count, :].todense()
        x += distance_weight[0] * sum([
            promoter_profiles[pid * sample_count:
                        (pid + 1) * sample_count, :].todense() for pid in pids
        ])
        x = x / float(len(promoters) + len(enhancers))
        xs.append(x)

    return xs

def gene_test(pg_pair, glist, ep_pair, sample_count, motif_count):
    rs = []
    ps = []
    cs = []
    ss = []
    for i, g in enumerate(glist):
        x = np.zeros((sample_count, motif_count))
        promoters = np.array(pg_pair[g])
        enhancers = [
            np.array(ep_pair[p]) for p in promoters if ep_pair.has_key(p)
        ]
        if len(enhancers) > 0:
            enhancers = np.unique(np.concatenate(enhancers))
        else:
            continue
        pids = np.array([
            np.where(promoter_results._labels == promoter)[0][0]
            for promoter in promoters
        ])
        for enhancer in enhancers:
            eid = np.where(enhancer_results._labels == enhancer)[0][0]
            weight = get_weight(
                get_distance(promoters[0], enhancer, ep_dis), distance_breaks,
                distance_weight)
            x += weight * enhancer_profiles[eid * sample_count:
                                      (eid + 1) * sample_count, :].todense()
        x += distance_weight[0] * sum([
            promoter_profiles[pid * sample_count:
                        (pid + 1) * sample_count, :].todense() for pid in pids
        ])
        x = x / float(len(promoters) + len(enhancers))
        print(i)
        r, p, s, c = test_fast(g, x, ncoefs=10)
        rs.append(x)
        ps.append(p)
        ss.append(s)
        cs.append(c)

    return rs, ps, ss, cs
