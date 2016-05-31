/*  plugins/missing2ref.c -- sets missing genotypes to reference allele.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include <getopt.h>
#include "bcftools.h"

bcf_hdr_t *in_hdr, *out_hdr;
int32_t *gts = NULL, mgts = 0;
int *arr = NULL, marr = 0;
uint64_t nchanged = 0;
int new_gt = bcf_gt_unphased(0);
int use_major = 0;

// binary mask array (number of elements = number of samples in the VCF)
//    0 - sample will not be processed
//    1 - sample will be processed
int *samples_to_process = NULL;

const char *about(void)
{
    return "Set missing genotypes (\"./.\") to ref or major allele (\"0/0\" or \"0|0\").\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Set missing genotypes\n"
        "Usage: bcftools +missing2ref [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -p, --phased     Set to \"0|0\" \n"
        "   -m, --major      Set to major allele \n"
        "   -s, --samples    Names of samples to process (or not, if prefixed with ^)\n"
        "\n"
        "Example:\n"
        "   bcftools +missing2ref in.vcf -- -p\n"
        "   bcftools +missing2ref in.vcf -- -p -m\n"
        "\n";
}


int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    int i;

    // process all samples in the VCF by default
    int process_all = 1;

    char *samples_str = NULL;
    int c;
    static struct option loptions[] =
    {
        {"phased",0,0,'p'},
        {"major",0,0,'m'},
        {"samples",required_argument,0,'s'},
        {0,0,0,0}
    };

    while ((c = getopt_long(argc, argv, "mps:?h", loptions, NULL)) >= 0)
    {
        switch (c) 
        {
            case 'p': new_gt = bcf_gt_phased(0); break;
            case 'm': use_major = 1; break;
            case 's':
                samples_str = optarg;
                process_all = 0;
                break;
            case 'h':
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }
    in_hdr  = in;
    out_hdr = out;

    samples_to_process = (int *) malloc (sizeof (int) * bcf_hdr_nsamples(in_hdr));

    if (process_all) {
        for (i = 0; i < bcf_hdr_nsamples(in_hdr); i++)
            samples_to_process[i] = 1;
    } else {
        int exclude = samples_str[0] == '^';

        // parse the comma-separated list of samples
        int nsamples;
        char **samples = hts_readlist(exclude ? &samples_str[1] : samples_str, 0, &nsamples);

        for (i = 0; i < nsamples; i++) {
            fprintf(stderr, "%s\n", samples[i]);
        }

        // check if all of the specified samples exist in the VCF header
        for (i = 0; i < nsamples; ++i) {
            if (bcf_hdr_id2int(in_hdr, BCF_DT_SAMPLE, samples[i]) == -1) {
                error("One of the samples is not present in the header: %s \n",
                      samples[i]);
            }
        }

        // initialize the sample mask
        for (i = 0; i < bcf_hdr_nsamples(in_hdr); i++) {
            samples_to_process[i] = exclude ? 1 : 0;
        }
        // update the mask for samples which will be processed
        for (i=0; i < nsamples; i++) {
            int pos = bcf_hdr_id2int(in_hdr, BCF_DT_SAMPLE, samples[i]);
            samples_to_process[pos] = exclude ? 0 : 1;
        }

        free(samples);
    }

    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    int ngts = bcf_get_genotypes(in_hdr, rec, &gts, &mgts);
    int i, changed = 0;
    
    // Calculating allele frequency for each allele and determining major allele
    // only do this if use_major is true
    int majorAllele = -1;
    int maxAC = -1;
    int an = 0;
    if(use_major == 1){
        hts_expand(int,rec->n_allele,marr,arr);
        int ret = bcf_calc_ac(in_hdr,rec,arr,BCF_UN_FMT);
        if(ret > 0){
            for(i=0; i < rec->n_allele; ++i){
                an += arr[i];
                if(*(arr+i) > maxAC){
                    maxAC = *(arr+i);
                    majorAllele = i;
                }
            }
        }
        else{
            fprintf(stderr,"Warning: Could not calculate allele count at position %d\n", rec->pos);
            exit(1);
        }

        // replacing new_gt by major allele
        if(bcf_gt_is_phased(new_gt))
            new_gt = bcf_gt_phased(majorAllele);
        else
            new_gt = bcf_gt_unphased(majorAllele);
    }

    // replace GTs of those samples that are supposed to be processed
    for (i = 0; i < ngts; i++) {
        if (gts[i] == bcf_gt_missing && samples_to_process[i / 2]) {
            gts[i] = new_gt;
            changed++;
        }
    }
    nchanged += changed;
    if ( changed ) bcf_update_genotypes(out_hdr, rec, gts, ngts);
    return rec;
}

void destroy(void)
{
    free(arr);
    fprintf(stderr,"Filled %"PRId64" REF alleles\n", nchanged);
    free(gts);
    free(samples_to_process);
}


