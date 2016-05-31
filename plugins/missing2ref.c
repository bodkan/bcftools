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
#include <htslib/khash_str2int.h>
#include <inttypes.h>
#include <getopt.h>
#include "bcftools.h"

bcf_hdr_t *in_hdr, *out_hdr;
int32_t *gts = NULL, mgts = 0;
int *arr = NULL, marr = 0;
uint64_t nchanged = 0;
int new_gt = bcf_gt_unphased(0);
int use_major = 0;

void *hdr_samples = NULL;
void *usr_samples = NULL;

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
        "   -p, --phased       Set to \"0|0\" \n"
        "   -m, --major        Set to major allele \n"
        "   -s, --samples      Names of samples to process \n"
        "\n"
        "Example:\n"
        "   bcftools +missing2ref in.vcf -- -p\n"
        "   bcftools +missing2ref in.vcf -- -p -m\n"
        "\n";
}


int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    int i;
    char *samples = NULL;
    int c;
    static struct option loptions[] =
    {
        {"phased",0,0,'p'},
        {"major",0,0,'m'},
        {"samples",required_argument,0,'s'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "mps:?h",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'p': new_gt = bcf_gt_phased(0); break;
            case 'm': use_major = 1; break;
            case 's': samples = optarg; break;
            case 'h':
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }
    in_hdr  = in;
    out_hdr = out;

    // load list of samples present in the VCF (specified in the header)
    hdr_samples = khash_str2int_init();
    for (i=0; i < bcf_hdr_nsamples(in_hdr); i++) {
        khash_str2int_inc(hdr_samples, bcf_hdr_int2id(in_hdr,BCF_DT_SAMPLE,i));
    }

    // parse the comma-separated list of samples for processing
    int user_nsamples;
    char **user_samples = hts_readlist(samples, 0, &user_nsamples);

    // check if alle the specified samples exist in the VCF header
    for (i = 0; i < user_nsamples; ++i) {
        if (!khash_str2int_has_key(hdr_samples, user_samples[i])) {
            error("One of the samples is not present in the header: %s \n",
                  user_samples[i]);
        }
    }

    // convert the list of user-specified samples to a hash table
    usr_samples = khash_str2int_init();
    for (i=0; i < user_nsamples; i++) {
        khash_str2int_inc(usr_samples, user_samples[i]);
    }

    fprintf(stderr, "Number of samples in the VCF: %d\n",  bcf_hdr_nsamples(in_hdr));
    fprintf(stderr, "Samples present in the VCF:\n");
    for (i = 0; i < bcf_hdr_nsamples(in_hdr); i++) {
        fprintf(stderr, "%5d %s\n", i, bcf_hdr_int2id(in_hdr, BCF_DT_SAMPLE, i));
    }

    fprintf(stderr, "\nNumber of samples specified by the user: %d\n", user_nsamples);
    fprintf(stderr, "Samples specified by the user:\n");
    for (i = 0; i < user_nsamples; i++) {
        fprintf(stderr, "%5d %s\n", i, user_samples[i]);
    }

    free(user_samples);

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

    // replace gts
    for (i=0; i<ngts; i++)
    {
        char *s = bcf_hdr_int2id(in_hdr, BCF_DT_SAMPLE, i / 2);
//        fprintf(stderr, "%d/%s/%d ", i / 2, s, khash_str2int_has_key(usr_samples, s));
//        fprintf(stderr, "%d ", gts[i]);
        if ( gts[i]==bcf_gt_missing && khash_str2int_has_key(usr_samples, s))
        {
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
}


