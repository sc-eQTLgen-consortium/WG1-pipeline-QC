#!/usr/bin/env python3
# Author: H.J.Westra (Edited by M.Vochteloo)
import argparse
import gzip
import os


parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-i", "--input", type=str, dest='input_path', required=True, help="The input VCF file.")
parser.add_argument("-o", "--output", type=str, dest='output_path', required=True, help="The output VCF file.")
parser.add_argument("-l", "--log", type=str, dest='log_path', required=True, help="The output logfile.")
parser.add_argument("-s", "--samples", nargs="*", type=str, dest='sample_filter_set', default=None, help="The samples to select.")
parser.add_argument("--sex", type=str, dest='sex_path', default=None, help="The sample-sex data file.")
parser.add_argument("-gq", "--genotype_quality", type=float, dest='thresh_gq', default=20., help="The genotype quality threshold. Default: 20.")
parser.add_argument("-abl", "--allelic_balance_lower", type=float, dest='thresh_ab_lower', default=0.2, help="The allelic balance lower-bound threshold. Default: 0.2.")
parser.add_argument("-abu", "--allelic_balance_upper", type=float, dest='thresh_ab_upper', default=0.8, help="The allelic balance upper-bound threshold. Default: 0.8.")
parser.add_argument("-ib", "--inbreeding_coefficient", type=float, dest='thresh_ib', default=-0.3, help="The inbreeding coefficient threshold. Default: -0.3.")
parser.add_argument("--no_snv_vqsr_check", dest='check_vqsr_snv', action='store_false', help="Add this flag to turn off the VQSR of SNV's. Default: True.")
parser.add_argument("--vqsr_snv", type=float, dest='thresh_vqsr_snv', default=99.8, help="The variant_quality_score_recalibration threshold. Default: 99.8.")
parser.add_argument("--no_indel_vqsr_check", dest='check_vqsr_indel', action='store_false', help="Add this flag to turn off the VQSR of indel's. Default: True.")
parser.add_argument("--vqsr_indel", type=float, dest='thresh_vqsr_indel', default=99.95, help="The variant_quality_score_recalibration threshold. Default: 99.95.")
parser.add_argument("--keep_multialleic", dest='remove_multiallelic', action='store_false', help="Add this flag to remove multiallelic variants. Default: True.")
parser.add_argument("--keep_non_pass_snv", dest='remove_non_pass_snv', action='store_false', help="Add this flag to remove SNV's without PASS label. Default: True.")
parser.add_argument("--keep_non_pass_indel", dest='remove_non_pass_indel', action='store_false', help="Add this flag to remove indel's without PASS label. Default: True.")
parser.add_argument("--keep_low_complexity", dest='filer_low_complexity', action='store_false', help="Add this flag to keep low complexity variants. Default: True.")
parser.add_argument("-maf", "--minor_allele_frequency", type=float, dest='thresh_maf', default=0.01, help="The minor allele frequency threshold. Default: 0.01.")
parser.add_argument("-cr", "--call_rate", type=float, dest='thresh_cr', default=0.99, help="The call rate threshold. Default: 0.99.")
parser.add_argument("-hwe", "--hardy_weinberg_equilibrium", type=float, dest='thresh_hwe', default=1E-6, help="The hardy weinberg equilibrium threshold. Default: 1e-6.")
parser.add_argument("-dp", "--filtered_depth", type=float, dest='thresh_dp', default=10., help="The filtered depth threshold. Default: 10.")
parser.add_argument("--keep_info_column", dest='strip_info_col', action='store_false', help="Add this flag to keep the INFO column. Default: True.")
args = parser.parse_args()

if not os.path.isdir(os.path.dirname(args.output_path)):
    os.mkdir(os.path.dirname(args.output_path))

sample_filter_set = args.sample_filter_set
if isinstance(sample_filter_set, list):
    sample_filter_set = set(sample_filter_set)

ismale_dict = None
if args.sex_path is not None:
    ismale_dict = {}
    fs = open(args.sex_path, 'r')
    iid_index = None
    sex_index = None
    for line_num, line in enumerate(fs):
        values = line.strip("\n").split("\t")
        if line_num == 0:
            for col_index, value in enumerate(values):
                if value == "IID":
                    iid_index = col_index
                elif value == "SEX":
                    sex_index = col_index
                else:
                    pass
            continue
        if values[sex_index] == "0" or values[sex_index] == "1":
            ismale_dict[values[iid_index]] = True
        elif values[sex_index] == "2":
            ismale_dict[values[iid_index]] = False
        else:
            print("Unexpected input in sex file.")
            exit()
    fs.close()


def parse_line(line_number, line, sample_mask, sample_male_mask):
    # line with genotype data
    elems = line.strip().split("\t", 10)

    # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    # 0     1       2       3       4       5       6       7       8
    chr = elems[0]
    pos = elems[1]

    id = elems[2]
    logid = chr + ":" + pos + ":" + id
    ref = elems[3]
    alt = elems[4].split(",")

    # check if this is a multiAllelic
    if args.remove_multiallelic and len(alt) > 1:
        return [line_number, False, logid + "\tMultiAllelic\t-\n"]

    # check if this variant is an indel
    is_indel = False
    if len(ref) > 1:
        is_indel = True

    if not is_indel:
        for a in alt:
            if len(a) > 1:
                is_indel = True

    filter = elems[6]
    # check VQSR first
    if is_indel:
        if filter.startswith("VQSR") and args.check_vqsr_indel:
            if not parse_vqsr(filter, is_indel):
                return [line_number, False, logid + "\tIndelBelowVQSR\t-\n"]
        elif args.remove_non_pass_indel and filter != "PASS":
            return [line_number, False, logid + "\tIndelNonPass\t-\n"]
    else:
        if filter.startswith("VQSR") and args.check_vqsr_snv:
            if not parse_vqsr(filter, is_indel):
                return [line_number, False, logid + "\tSNVBelowVQSR\t-\n"]
        elif args.remove_non_pass_snv and filter != "PASS":
            return [line_number, False, logid + "\tSNVNonPass\t-\n"]

    # parse info string
    infoelems = elems[7].split(";")
    for elem in infoelems:
        vals = elem.split("=")
        if len(vals) > 1:
            key = vals[0]
            value = vals[1]
            if key == "InbreedingCoeff":
                if value == ".":
                    return [line_number, False, logid + "\tIncorrectInbreedingCoeff:" + value + "\t-\n"]
                valuef = float(value)
                if valuef < args.thresh_ib:
                    return [line_number, False, logid + "\tBelowInbreedingCoeff:" + value + "\t-\n"]
        # other filters based on info column can be implemented here...

    # parse genotypes
    elems = line.strip().split("\t")
    format = elems[8].split(":")
    gtcol = -1
    adcol = -1
    dpcol = -1
    gqcol = -1
    for c in range(0, len(format)):
        if format[c] == "GT":
            gtcol = c
        elif format[c] == "AD":
            adcol = c
        elif format[c] == "DP":
            dpcol = c
        elif format[c] == "GQ":
            gqcol = c

    format_out = ""
    if gtcol > -1:
        format_out = "GT"
    if dpcol > -1:
        format_out += ":DP"
    if adcol > -1:
        format_out += ":AD"
        format_out += ":AB"
    if gqcol > -1:
        format_out += ":GQ"

    # only continue if there are genotypes to parse
    if gtcol < 0:
        return [line_number, False, logid + "\tNoGTCol\t-\n"]

    # get current dosages
    # store split columns temporarily
    sampledataelems = [None] * len(elems)
    dosages_pre_filter = []
    dosage_pre_filter = None
    sampledosages = [-1] * len(elems)
    cctr = 0
    gtsep = ""
    for i in range(9, len(elems)):
        if sample_mask[cctr]:
            # split the sample columns
            sampledata = elems[i].split(":")
            sampledataelems[i] = sampledata
            gt = sampledata[gtcol].split("/")
            if len(gt) < 2:
                # less than two alleles; perhaps the data is phased
                gt = sampledata[gtcol].split("|")
                gtsep = "|"
            if len(gt) < 2 or len(gt) > 2:
                # malformatted genotypes
                # more than two alleles.. skip variant?
                # fhlog.write(headerInfo[1]+"\tPolyPloid?\t-\n"
                dosages_pre_filter.append(-1)
            else:
                if gt[0] == ".":
                    dosages_pre_filter.append(-1)
                else:
                    dosage_pre_filter = int(gt[0]) + int(gt[1])
                    dosages_pre_filter.append(dosage_pre_filter)
                    sampledosages[i] = dosage_pre_filter
        cctr += 1

    # calculate allele frequency etc before filtering
    stats = get_stats(
        dosages=dosages_pre_filter,
        sample_male_mask=sample_male_mask,
        chr=chr
    )
    if not stats[0]:
        return [line_number, False, logid + "\tFailedPrefilterVarStats\t" + stats[1] + "\n"]

    # check whether variant becomes monomorphic after filtering poor calls
    # iterate samples
    dosages_post_filter = []
    sampleinfos = []
    nr_genotypes_replaced = 0
    average_depth = 0
    average_depth_calls = 0
    poor_dp = 0
    poor_gq = 0
    poor_ab_hom_a = 0
    poor_ab_hom_b = 0
    poor_ab_het = 0
    cctr = 0
    for i in range(9, len(elems)):
        sampledata = sampledataelems[i]
        if sampledata is not None:
            dosage = sampledosages[i]
            ab = -1
            dp = 0
            # if the genotype is not missing
            # check genotype qual
            if dpcol > -1 and dpcol < len(sampledata):  # ugh, VCF really allows all kind of stuff in their genotype format... even if you parse the FORMAT column, there's no telling what is actually provided per sample..
                # determine read depth
                if sampledata[dpcol] == ".":
                    dp = 0
                else:
                    dp = int(sampledata[dpcol])
                    average_depth += dp
                    average_depth_calls += 1
                if dp < args.thresh_dp:
                    dosage = -1
                    poor_dp += 1
                    nr_genotypes_replaced += 1
            if gqcol > -1 and dosage > -1 and gqcol < len(sampledata):
                gq = float(sampledata[gqcol])
                if gq < args.thresh_gq:
                    if dosage == 0 and gq == 0 and dp >= args.thresh_dp:
                        # bypass potential error in gVCF merging, see: https://github.com/broadinstitute/gatk/issues/5445
                        # best way would to also parse the PL field, if available, and check whether the homRef and het fields are both 0 (unlikely for a homRef call)
                        pass
                    else:
                        dosage = -1
                        nr_genotypes_replaced += 1
                        poor_gq += 1
            # check allelic balance
            # only need to check this for hets, I guess
            if dosage > -1 and adcol > -1 and adcol < len(sampledata):
                ad = sampledata[adcol].split(",")
                ad1 = float(ad[0])
                ad2 = float(ad[1])
                if ad1 + ad2 == 0:
                    dosage = -1
                    nr_genotypes_replaced += 1
                else:
                    ab = ad2 / (ad1 + ad2)
                    # data is already corrected for DP, so no need to check again.
                    if dosage == 0:
                        if ab > args.thresh_ab_lower:
                            dosage = -1
                            nr_genotypes_replaced += 1
                            poor_ab_hom_a += 1
                    if dosage == 1:
                        if ab < args.thresh_ab_lower or ab > args.thresh_ab_upper:
                            dosage = -1
                            nr_genotypes_replaced += 1
                            poor_ab_het += 1
                    if dosage == 2:
                        if ab < args.thresh_ab_upper:
                            dosage = -1
                            nr_genotypes_replaced += 1
                            poor_ab_hom_b += 1
            dosages_post_filter.append(dosage)

            sampleinfo = []
            gtout = sampledata[gtcol]
            if dosage_pre_filter == -1:
                gtout = "." + gtsep + "."
            sampleinfo.append(gtout)

            # append original information, but note that missing genotypes
            # do not always conform to the FORMAT string
            if dpcol > -1 and dpcol < len(sampledata):
                sampleinfo.append(sampledata[dpcol])
            if adcol > -1 and adcol < len(sampledata):
                sampleinfo.append(sampledata[adcol])
                sampleinfo.append(str(round(ab, 2)))
            if gqcol > -1 and gqcol < len(sampledata):
                sampleinfo.append(sampledata[gqcol])

            sampleinfos.append(":".join(sampleinfo))
        cctr += 1

    # average read depth
    if average_depth_calls > 0:
        average_depth /= average_depth_calls
    else:
        average_depth = 0

    # variant level quals
    _, _, _, _, _, maf_postfilter = get_maf(
        dosages=dosages_post_filter,
        sample_male_mask=sample_male_mask,
        chr=chr
    )
    if maf_postfilter == 0:
        logoutln = "{}\tMonomorphicPostFilter\t{};NrGenosReplaced:{}" \
                   ";PoorDP:{};AvgDP:{:.2f} ({} calls);PoorGQ:{}" \
                   ";PoorABHomA:{};PoorABHomB:{};PoorABHet:{}" \
                   "\n".format(logid,
                               stats[1],
                               nr_genotypes_replaced,
                               poor_dp,
                               average_depth,
                               average_depth_calls,
                               poor_gq,
                               poor_ab_hom_a,
                               poor_ab_hom_b,
                               poor_ab_het)
        return [line_number, False, logoutln]

    format_out = ""
    if gtcol > -1:
        format_out = "GT"
    if dpcol > -1:
        format_out += ":DP"
    if adcol > -1:
        format_out += ":AD"
        format_out += ":AB"
    if gqcol > -1:
        format_out += ":GQ"

    elems[8] = format_out

    if args.strip_info_col:
        elems[7] = "."
    outln = "\t".join(elems[0:9])  # construct line header
    # write variant, somehow
    outln += "\t" + "\t".join(sampleinfos) + "\n"
    return [line_number, True, logid + "\tPASSQC\t" + stats[1] + "\n", outln]


def parse_vqsr(vqsrstring, is_indel):
    if not vqsrstring.startswith("VQSR"):
        return True
    vqsrstring = vqsrstring.replace("VQSRTrancheINDEL", "")
    vqsrstring = vqsrstring.replace("VQSRTrancheSNP", "")
    vqsrstring = vqsrstring.replace("+", "")
    velems = vqsrstring.split("to")
    lower = float(velems[0])
    if is_indel:
        return lower >= args.thresh_vqsr_indel
    else:
        return lower >= args.thresh_vqsr_snv


def get_stats(dosages, sample_male_mask, chr):
    nrcalled, nrdosages, nrhoma, nrhets, nrhomb, maf = get_maf(
        dosages=dosages,
        sample_male_mask=sample_male_mask,
        chr=chr
    )
    if nrcalled == 0:
        return [False, "CR=0.000;MAF=0.000;HWE=-1.000e+00"]

    var_ok = True
    callrate = nrcalled / nrdosages
    hwe = -1
    if callrate >= args.thresh_cr and maf >= args.thresh_maf:
        hwe = calculate_hwe(
            obs_hets=nrhets,
            obs_hom1=nrhoma,
            obs_hom2=nrhomb
        )  # only calculate HWE when MAF and CR above thresholds
        if hwe < args.thresh_hwe:
            var_ok = False
    else:
        var_ok = False
    statout = "CR={:.3f};MAF={:.3f};HWE={:0.3e}".format(callrate, maf, hwe)

    return [var_ok, statout]


def get_maf(dosages, sample_male_mask, chr):
    nrhoma = 0
    nrhets = 0
    nrhomb = 0
    nrmiss = 0

    for dosage, is_male in zip(dosages, sample_male_mask):
        if chr in ["X", "y"] and (is_male is None or is_male):
            continue

        if dosage == 0:
            nrhoma += 1
        elif dosage == 1:
            nrhets += 1
        elif dosage == 2:
            nrhomb += 1
        else:
            nrmiss += 1

    nrcalled = nrhoma + nrhets + nrhomb
    nrdosages = nrcalled + nrmiss
    maf = 0
    if nrcalled > 0:
        maf = (nrhoma * 2 + nrhets) / (nrcalled * 2)
        if maf > 0.5:
            maf = 1 - maf

    return nrcalled, nrdosages, nrhoma, nrhets, nrhomb, maf


def calculate_hwe(obs_hets, obs_hom1, obs_hom2):
    obs_homc = obs_hom1
    if obs_hom1 < obs_hom2:
        obs_homc = obs_hom2
    obs_homr = obs_hom2
    if obs_hom1 < obs_hom2:
        obs_homr = obs_hom1

    rare_copies = 2 * obs_homr + obs_hets
    l_genotypes = obs_hets + obs_homc + obs_homr

    if l_genotypes == 0:
        return -1

    het_probs = [0] * (rare_copies + 1)
    mid = int(rare_copies) * int(2 * l_genotypes - rare_copies) / int(2 * l_genotypes)
    mid = int(mid)
    if mid % 2 != rare_copies % 2:
        mid += 1
    mid = int(mid)

    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = l_genotypes - curr_hets - curr_homr
    het_probs[int(mid)] = 1.0
    sum = het_probs[int(mid)]

    curr_hets = int(mid)
    while curr_hets > 1:
        het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        sum += het_probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1
        curr_hets -= 2

    curr_hets = int(mid)
    curr_homr = (rare_copies - mid) / 2
    curr_homc = l_genotypes - curr_hets - curr_homr
    while curr_hets <= (rare_copies - 2):
        het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        sum += het_probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1
        curr_hets += 2

    i = 0
    while i <= rare_copies:
        het_probs[i] /= sum
        i += 1

    p_hwe = 0.0
    i = 0
    while i <= rare_copies:
        if het_probs[i] <= het_probs[obs_hets]:
            p_hwe += het_probs[i]
        i += 1

    if p_hwe > 1:
        p_hwe = 1
    return p_hwe

fh = gzip.open(args.input_path, 'rt')
fho = gzip.open(args.output_path, 'wt')
fhlog = gzip.open(args.log_path, 'wt')
fhlog.write("Id\tReason\tStats\n")

sample_mask = []
sample_male_mask = []

linectr = 0
lineswritten = 0
for line in fh:
    if line.startswith("#CHROM"):
        outln = "##{} tresh_GQ={};tresh_AB_lower={};" \
                "tresh_AB_upper={};tresh_IB={};check_VQSR_SNV={};" \
                "tresh_VQSR_SNV={};check_VQSR_Indel={};" \
                "tresh_VQSR_Indel={};remove_multiAllelic={};" \
                "remove_nonPASS_SNV={};remove_nonPASS_indel={};" \
                "filter_lowcomplexity={};tresh_MAF={};tresh_CR={};" \
                "tresh_HWE={};tresh_DP={}" \
                "strip_INFO_col={}\n".format(args.input_path,
                                             args.thresh_gq,
                                             args.thresh_ab_lower,
                                             args.thresh_ab_upper,
                                             args.thresh_ib,
                                             args.check_vqsr_snv,
                                             args.thresh_vqsr_snv,
                                             args.check_vqsr_indel,
                                             args.thresh_vqsr_indel,
                                             args.remove_multiallelic,
                                             args.remove_non_pass_snv,
                                             args.remove_non_pass_indel,
                                             args.filer_low_complexity,
                                             args.thresh_maf,
                                             args.thresh_cr,
                                             args.thresh_hwe,
                                             args.thresh_dp,
                                             args.strip_info_col
                                             )
        fho.write(outln)

        # header line with samples
        elems = line.strip().split("\t")
        elems[6] = "FILTER"
        elems[7] = "INFO"
        headerout = "\t".join(elems[0:9])
        n_samples_included = 0
        n_samples_excluded = 0
        for i in range(9, len(elems)):
            sample = elems[i]
            if sample_filter_set is None or sample in sample_filter_set:
                headerout += "\t" + sample
                sample_mask.append(True)
                sample_male_mask.append(None if ismale_dict is None else ismale_dict[sample])
                n_samples_included += 1
            else:
                sample_mask.append(False)
                n_samples_excluded += 1
        print("{:,} samples selected from VCF header, {:,} excluded".format(n_samples_included, n_samples_excluded))
        fho.write(headerout + "\n")
        lineswritten += 1
    elif line.startswith("#"):
        if line.startswith("##FORMAT") and line[13:15] not in ["GT", "DP", "AD", "AB", "GQ"]:
            continue
        if line.startswith("##INFO") and args.strip_info_col:
            continue
        fho.write(line)
        lineswritten += 1
    else:
        # this can now be potentially multithreaded?
        parsed = parse_line(
            line_number=linectr,
            line=line,
            sample_mask=sample_mask,
            sample_male_mask=sample_male_mask
        )  # returns: [line_number, False, logid+"\PASSQC\t"+stats[1], outln]
        fhlog.write(parsed[2])
        if parsed[1]:
            fho.write(parsed[3])
            lineswritten += 1
    linectr += 1
    if linectr % 10000 == 0:
        print("{:,} lines parsed, {:,} written".format(linectr, lineswritten), end='\r')
        fho.flush()
        fhlog.flush()
print("{:,} lines parsed, {:,} written".format(linectr, lineswritten), end='\n')
print("Done. How about that!")
fh.close()
fho.close()
fhlog.close()