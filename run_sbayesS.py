import hail as hl
import argparse
import hailtop.batch as hb
from hailtop.batch.job import BashJob


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
    
    sumstats_path = f'gs://gwasns-analysis/panukbb/sumstats/{trait}.ma'
    return trait, sumstats_path


def run_sbayesS(b: hb.batch.Batch,
                image: str,
                seq: int
                ):
    j = b.new_job(name=f'run-sbayesS-{seq}')

    j.image(image)
    j.cpu(4)
    j.memory('highmem')

    # input for sbayesS
    trait = format_to_ma(seq)[0]
    sumstats_path = format_to_ma(seq)[1]
    ldref_path = f'gs://gwasns-analysis/panukbb/ldref/ukbEURu_hm3_sparse_mldm_list.txt'
    out_path = f'gs://gwasns-analysis/panukbb/outputs/{trait}'

    # not sure about this part?
    j.declare_resource_group(outf = {
        'mcmcsamples.Par' : '{root/sbayesS}.mcmcsamples.Par',
        'mcmcsamples.SnpEffects' : '{root/sbayesS}.mcmcsamples.SnpEffects',
        'parRes' : '{root/sbayesS}.parRes',
        'snpRes' : '{root/sbayesS}.snpRes'
    })

    j.command(f'''
    gctb --sbayes S --mldm  {ldref_path} \
    --gwas-summary {sumstats_path} \
    --chain-length 10000 --burn-in 2000 --out-freq 10 \
    --robust --exclude-mhc --num-chains 4 \
    --out {j.outf}''')

    b.write_output(j.outf, f'{out_path}')

    return j


def main(args):
    backend = hb.ServiceBackend(billing_project='ukbb_diverse_pops',
                                bucket='ukbb-diverse-pops')
    b = hb.Batch(backend=backend, name='sbayesS')
    sbayesS_img = 'gcr.io/ukbb-diversepops-neale/ywang-sbrs:test'

    sbayesS_jobs = []
    for idx in range(args.idx_start, args.idx_end, args.idx_step):
        sbayesS_j = run_sbayesS(b = b, image = sbayesS_img, seq = idx)
        sbayesS_jobs.append(sbayesS_j)

    b.run(open = True)
    backend.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--idx_start', type = int, required=True, help='The start index for QCd phenotype, ranged from 1:520 in EUR')
    parser.add_argument('--idx_end', type = int, required=True,
                        help='The end index for QCd phenotype, ranged from 1:520 in EUR')
    parser.add_argument('--idx_step', type = int, required=True,
                        help='The step index for QCd phenotype, ranged from 1:520 in EUR')
    args = parser.parse_args()

    main(args)
