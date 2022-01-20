import hailtop.batch as hb
import argparse

def run_format(b: hb.batch.Batch,
               image: str,
               seq: int):
    cores = 4
    g = b.new_job(name='formatGWAS')
    g.image(image)
    g.cpu(cores)
    g.command(f'''
    python3 /sbayesS/format_to_ma.py \
        --seq {seq} 
    ''')
    return g


def run_sbayesS(b: hb.batch.Batch,
                image: str,
                trait: str,
                depends_on_j
                ):
    j = b.new_job(name=f'run-sbayesS-{trait}')
    j.depends_on(depends_on_j)
    j.image(image)
    j.cpu(4)
    j.memory('highmem')

    # input for sbayesS
    sumstats_path = f'gs://gwasns-analysis/panukbb/sumstats/{trait}.ma'
    ldref_path = f'gs://gwasns-analysis/panukbb/ldref/ukbEURu_hm3_sparse_mldm_list.txt'
    out_path = f'gs://gwasns-analysis/panukbb/outputs/{trait}'

    # not sure about this part?
    j.declare_resource_group(outf={
        'mcmcsamples.Par': '{root}.mcmcsamples.Par',
        'mcmcsamples.SnpEffects': '{root}.mcmcsamples.SnpEffects',
        'parRes': '{root}.parRes',
        'snpRes': '{root}.snpRes'
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

    format_jobs = []
    sbayesS_jobs = []
    for idx in range(args.idx_start, args.idx_end, args.idx_step):
        format_g = run_format(b=b, image=sbayesS_img, seq=idx)
        format_jobs.append(format_g)

        sbayesS_j = run_sbayesS(b=b, image=sbayesS_img, depends_on_j=format_g, trait=format_g.trait)
        sbayesS_jobs.append(sbayesS_j)

    b.run(open=True)
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
