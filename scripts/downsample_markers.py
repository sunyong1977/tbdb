import sys
import csv
import argparse
from collections import defaultdict
import json


def main(args):
    lineage_markers = defaultdict(set)
    lineages = set()
    for row in csv.DictReader(open(args.csv)):
        lineages.add(row["Lin"])
        if len(lineage_markers[row["Lin"]])>=10:
            continue
        if row["essential"]=="1" and row["dr_bin"]=="0" and row["pe_ppe_bin"]=="0" and row["syn_nonsyn"]=="synonymous":
            lineage_markers[row["Lin"]].add(json.dumps(row))


    for row in csv.DictReader(open(args.csv)):
        if len(lineage_markers[row["Lin"]])>=10:
            continue
        if row["dr_bin"]=="0" and row["pe_ppe_bin"]=="0" and row["syn_nonsyn"]=="synonymous":
            lineage_markers[row["Lin"]].add(json.dumps(row))

    for row in csv.DictReader(open(args.csv)):
        if len(lineage_markers[row["Lin"]])>=10:
            continue
        if row["dr_bin"]=="0" and row["pe_ppe_bin"]=="0":
            lineage_markers[row["Lin"]].add(json.dumps(row))

    final_markers = []
    for lin in lineages:
        for marker in lineage_markers[lin]:
            tmp = json.loads(marker)
            tmp["Pos"] = int(tmp["Pos"])
            final_markers.append(tmp)

    with open(args.out,"w") as O:
        for marker in sorted(final_markers,key=lambda x: x["Pos"]):
            allele = marker["Alt"] if (marker["Lin"]!="4" and marker["Lin"]!="4.9") else marker["Ref"]
            O.write("Chromosome\t%s\t%s\tlineage%s\t%s\n" % (marker["Pos"]-1, marker["Pos"], marker["Lin"],allele))


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--csv',type=str,help='File with samples')
parser.add_argument('--out',type=str,help='File with samples')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
