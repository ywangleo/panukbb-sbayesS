import hail as hl
import hailtop.batch as hb
import argparse


def format_to_ma(trait: str,
                 Ns: str):
    # HM3 SNP list for SBayesS
    # snps = hl.import_table("gs://gwasns-analysis/panukbb/ukbEURu_hm3_v3_50k.ldm.sparse.info")
    # snps = snps.annotate(locus=snps.Chrom + ":" + snps.PhysPos)
    # snps = snps.key_by("locus")
    #
    # # format GWAS
    # trait_path = f'gs://ukb-diverse-pops/ld_prune/export_results/2112_final_results/{trait}.tsv.bgz'
    # gwas = hl.import_table(trait_path)
    # gwas = gwas.annotate(locus=(gwas.chr + ":" + gwas.pos), N=Ns)
    # gwas = gwas.key_by("locus")
    # gwas = gwas.semi_join(snps)
    # gwas2 = gwas.annotate(SNP=snps[gwas.locus].SNP)
    #
    # gwas3 = gwas2.filter(hl.is_defined(gwas2.SNP), keep=True)
    # gwas3 = gwas3.key_by()
    # gwas4 = gwas3.select(gwas3.SNP, gwas3.ref, gwas3.alt, gwas3.af_EUR, gwas3.beta_EUR, gwas3.se_EUR, gwas3.pval_EUR,
    #                     gwas3.N)
    # gwas4 = gwas4.rename(
    #     {'ref': 'A1', 'alt': 'A2', 'af_EUR': 'freq', 'beta_EUR': 'b', 'se_EUR': 'se', 'pval_EUR': 'p'})
    # gwas4.export(f'gs://gwasns-analysis/panukbb/sumstats/{trait}.ma')
    return trait


def run_sbayesS(b: hb.batch.Batch,
                image: str,
                trait: str,
                depends_on_j
                ):
    j = b.new_job(name=f'run-sbayesS-{trait}')
    j.depends_on(depends_on_j)
    j.image(image)
    j.cpu(1)
    j.memory('highmem')
    j.storage(f'45Gi')

    # input for sbayesS
    sumstats_path = b.read_input(f'gs://gwasns-analysis/panukbb/sumstats/{trait}.ma')
    out_path = f'gs://gwasns-analysis/panukbb/outputs/{trait}'
    ldref_dir = b.read_input('gs://gwasns-analysis/panukbb/ldref/')

    # outputs for sbayesS
    j.declare_resource_group(outf={
        'mcmcsamples.Par': '{root}.mcmcsamples.Par',
        'mcmcsamples.SnpEffects': '{root}.mcmcsamples.SnpEffects',
        'parRes': '{root}.parRes',
        'snpRes': '{root}.snpRes'
    })

    j.command(f'''
    python3.8 -c '
with open("io/ukbEURu_hm3_sparse_mldm_list.txt.remapped", "w") as new_mldm_list:
    with open("{ldref_dir}/ukbEURu_hm3_sparse_mldm_list.txt") as old_mldm_list:
        orig_file_names = [line.rstrip() for line in old_mldm_list]
        for orig_path in orig_file_names:
            new_mldm_list.write(orig_path.replace("gs://gwasns-analysis/panukbb/ldref", """{ldref_dir}""") + """\n""")
    '
    
    /sbayesS/gctb_2.03beta_Linux/gctb --sbayes S --mldm /io/ukbEURu_hm3_sparse_mldm_list.txt.remapped \
        --gwas-summary {sumstats_path} \
        --chain-length 10000 --burn-in 2000 --out-freq 10 \
        --robust --exclude-mhc --num-chains 4 \
        --out {j.outf}

    ''')

    b.write_output(j.outf, f'{out_path}')


def main(args):
    backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                remote_tmpdir='gs://ukb-diverse-pops')

    b = hb.Batch(backend=backend, name='sbayesS')
    sbayesS_img = 'gcr.io/ukbb-diversepops-neale/ywang-sbrs:test'

    #image = hb.build_python_image('gcr.io/ukbb-diversepops-neale/batch-python:test',
    #                           requirements=['hail'])

    # extract sample size information
    phenos = hl.import_table('gs://gwasns-analysis/panukbb/panukbbEUR_qcd_phenos_Ns.txt')

    for idx in range(args.idx_start, args.idx_end, args.idx_step):
        tmp = phenos.filter(hl.int(phenos.idx) == idx)
        trait = tmp.phenotype_id.collect()[0]
        Ns = tmp.Ntotal.collect()[0]

        format_job = b.new_python_job(name=f'format-{trait}')
        format_job.image('gcr.io/ukbb-diversepops-neale/batch-python:test').call(format_to_ma, trait=trait, Ns=Ns)
        run_sbayesS(b=b, image=sbayesS_img, depends_on_j=format_job, trait=trait)

    b.run()
    backend.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--idx_start', type=int, required=True,
                        help='The start index for QCd phenotype, ranged from 1:520 in EUR')
    parser.add_argument('--idx_end', type=int, required=True,
                        help='The end index for QCd phenotype, ranged from 1:520 in EUR')
    parser.add_argument('--idx_step', type=int, required=True,
                        help='The step index for QCd phenotype, ranged from 1:520 in EUR')
    args = parser.parse_args()

    main(args)
