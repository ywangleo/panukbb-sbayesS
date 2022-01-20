# pip install ukbb-common==0.1.2
# git clone https://github.com/atgu/ukbb_pan_ancestry
# git clone https://github.com/Nealelab/ukb_common
# curl -sSL https://broad.io/install-gcs-connector | python
from ukbb_pan_ancestry.resources.results import get_variant_results_path

from ukbb_pan_ancestry.export_results import get_pheno_id
from ukbb_pan_ancestry.heritability.import_heritability import get_h2_ht
import hail as hl

hl.init(spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',
                    'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})
from ukbb_pan_ancestry import *

# get the full QCd phenotype list (a total of 520 phenotypes)
h2_qc_ht = hl.read_table(get_h2_ht())
h2_qc_ht = h2_qc_ht.explode('heritability')
h2_qc_ht1 = h2_qc_ht.filter(h2_qc_ht.heritability.pop == 'EUR')
h2_qc_ht2 = h2_qc_ht1.filter((h2_qc_ht1.heritability.qcflags.pass_all == True))  ##520 phenotypes
h2_qc_ht2 = h2_qc_ht2.annotate(phenotype_id=get_pheno_id(tb=h2_qc_ht2))
h2_qc_ht2 = h2_qc_ht2.key_by(h2_qc_ht2.phenotype_id)

# extract sample size information for EUR based on QCd phenotype list
pheno_manifest = hl.read_matrix_table(get_variant_results_path('full')).cols()
pheno_manifest2 = pheno_manifest.annotate(phenotype_id=get_pheno_id(tb=pheno_manifest))
pheno_manifest2 = pheno_manifest2.key_by(pheno_manifest2.phenotype_id)
pheno_manifest3 = pheno_manifest2.semi_join(h2_qc_ht2)  ##rows in h2_qc_ht2 based on matched key
pheno_manifest3 = pheno_manifest3.explode('pheno_data')
pheno_manifest3 = pheno_manifest3.filter(pheno_manifest3.pheno_data.pop == "EUR")


def get_Ns(ncase, ncontrol):
    return hl.if_else(hl.is_defined(ncontrol),
                      ncase + ncontrol,
                      ncase)


pheno_manifest3 = pheno_manifest3.annotate(
    Ntotal=get_Ns(pheno_manifest3.pheno_data.n_cases, pheno_manifest3.pheno_data.n_controls))
pheno_manifest3 = pheno_manifest3.annotate(prev=pheno_manifest3.pheno_data.n_cases / pheno_manifest3.Ntotal)

pheno_manifest4 = pheno_manifest3.key_by()
pheno_manifest4 = pheno_manifest4.select(pheno_manifest4.phenotype_id, pheno_manifest4.Ntotal, pheno_manifest4.prev)
pheno_manifest4 = pheno_manifest4.add_index()
pheno_manifest4.export("gs://gwasns-analysis/panukbb/panukbbEUR_qcd_phenos_Ns.txt")
