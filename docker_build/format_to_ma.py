import hail as hl
import argparse

def format_to_ma(seq: int):
    # extract sample size information
    phenos = hl.import_table('gs://gwasns-analysis/panukbb/panukbbEUR_qcd_phenos_Ns.txt')
    tmp = phenos.filter(hl.int(phenos.idx) == seq)
    trait = tmp.phenotype_id.collect()[0]
    Ns = tmp.Ntotal.collect()[0]

    # HM3 SNP list for SBayesS
    snps = hl.import_table("gs://gwasns-analysis/panukbb/ldref/ukbEURu_hm3_v3_50k.ldm.sparse.info")
    snps = snps.annotate(locus=snps.Chrom + ":" + snps.PhysPos)
    snps = snps.key_by("locus")

    # format GWAS
    trait_path = f'gs://ukb-diverse-pops/ld_prune/export_results/2112_final_results/{trait}.tsv.bgz'
    gwas = hl.import_table(trait_path)
    gwas = gwas.annotate(locus=(gwas.chr + ":" + gwas.pos), N=Ns)
    gwas = gwas.key_by("locus")
    gwas = gwas.semi_join(snps)
    gwas2 = gwas.annotate(SNP=snps[gwas.locus].SNP)

    gwas3 = gwas2.filter(hl.is_defined(gwas2.SNP), keep=True)
    gwas3 = gwas3.key_by()
    gwas4 = gwas3.select(gwas3.SNP, gwas3.ref, gwas3.alt, gwas3.af_EUR, gwas3.beta_EUR, gwas3.se_EUR, gwas3.pval_EUR,
                         gwas3.N)
    gwas4 = gwas4.rename(
        {'ref': 'A1', 'alt': 'A2', 'af_EUR': 'freq', 'beta_EUR': 'b', 'se_EUR': 'se', 'pval_EUR': 'p'})
    gwas4.export(f'gs://gwasns-analysis/panukbb/sumstats/{trait}.ma')

    return trait

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pheno_idx', type = int, required=True)
    args = parser.parse_args()

    format_to_ma(args.pheno_idx)
