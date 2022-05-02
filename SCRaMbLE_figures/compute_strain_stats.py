#!/usr/bin/env python

import argparse
import csv
from collections import defaultdict

DESCRIPTION_FILE="desc/create_adjacency_list.txt"

def load_description():
    lines = open(DESCRIPTION_FILE).readlines()
    return "\n".join(lines)

def get_halfsites(num_segments):
    hsites = []
    
    for i in range(1,num_segments+1):
            prefix = ""
            if abs(i) < 10:
                prefix = "0"
        
            hsites.append("%s%dL" % (prefix,i))
            hsites.append("%s%dR" % (prefix,i))
    
    return hsites
    
def create_junctions(num_segments):
    junctions = defaultdict(list)
    
    for i in range(1,num_segments+1):
        prefix = ""
        if abs(i) < 10:
            prefix = "0"
        
        junctions[i].append("%s%dL" % (prefix,i))
        junctions[i].append("%s%dR" % (prefix,i))
        
    print len(junctions.keys())
    return junctions    


def create_split_reconstruction(filename, junctions):
    rec = defaultdict(list)
    for record in open(filename,'r'):
        #useful to add comments
        if record[0] == '#':
            continue

        fields = record.strip().split("\t")
        nodes = fields[1].split(",")

        #calculate copy number
        copy_num = defaultdict(int)
        for node in nodes:
            copy_num[abs(int(node))] += 1
        
        for i in range(0, len(nodes)):
            segment = int(nodes[i])
            ends = junctions[abs(segment)]
            
            if segment < 0:
                rec[fields[0]].append((ends[1], ends[0], copy_num[abs(segment)]))
            else:
                rec[fields[0]].append((ends[0], ends[1], copy_num[abs(segment)]))
    
    return rec
    

def subset_segment(hsites,hstart, hstop):
    subset1 = set()
    s1 = int(hstart[:-1])
    s2 = int(hstop[:-1])
    
    start = min(s1,s2)
    stop  = max(s1,s2)
    
    return set(range(s1+1, s2))
   

def get_junction_type(h1,h2,c1,c2):
    hsites = get_halfsites(43)
    essential_segments = set([2,7,9,10,12,20])
    
    seg1 = h1
    seg2 = h2
    seg1_cn = c1
    seg2_cn = c2
    s1 = int(seg1[:-1])
    d1 = seg1[-1]
    s2 = int(seg2[:-1])
    d2 = seg2[-1]

    junction_type = None
    ss = subset_segment(hsites, seg1, seg2)
    
    if seg1_cn != seg2_cn:
        junction_type = 'duplicated junction'
    elif  (seg1 == "43R" and seg2 == "01L") or (seg1 == "01L" and seg2 == "43R"):
        junction_type = "parental"
    elif s2 == s1:
        junction_type = "repeat internal junction"
    elif abs(s2 - s1) == 1 and d1 != d2:
        junction_type = "parental"
    elif d1 == d2:
        junction_type = "inversion"    
    elif abs(s2 - s1) > 1:
        if ss.intersection(essential_segments):
            junction_type = "complex"
        else:
            junction_type = "deletion"
    else:
        print "unkown", seg1, seg2

    #print '{}({}), {}({}): {}'.format(seg1, seg1_cn, seg2, seg2_cn, junction_type)
    return junction_type    


def compute_halfsites_stats(rec):
    stats = defaultdict(lambda: defaultdict(int))
    
    for strain, genome in rec.iteritems():
        unique_junctions = set()
        for i in range(len(genome)):
            seg1 = genome[i][1]
            seg2 = genome[(i+1) % len(genome)][0]
            seg1_cn = genome[i][2]
            seg2_cn = genome[(i+1) % len(genome)][2]
            
            if ((seg1, seg2) not in unique_junctions) and ((seg2, seg1) not in unique_junctions):
                junction_type = get_junction_type(seg1,seg2,seg1_cn,seg2_cn)
                stats[strain][junction_type] += 1
                unique_junctions.add((seg1,seg2))
    return stats
            
def save_to_file(filename, stats):
    
    fh = csv.writer(open(filename, 'w'))
    fields = ['parental', 'deletion', 'inversion', 'complex', 'duplicated junction', 'repeat internal junction']
    header = ['strain']
    header.extend(fields)
    fh.writerow(header)
    hsites = get_halfsites(43)
    for h1 in stats.keys():
        record = [h1]
        for f in fields:
            record.append(str(stats[h1][f]))
        fh.writerow(record)
        
    
def parse_options():
    parser = argparse.ArgumentParser(description=load_description())
    parser.add_argument("-i", "--input-file", dest="input_filename", help="Reconstruction filename", type=str, required=True)
    parser.add_argument("-o", "--output-file", dest="output_filename", help="Output filename", type=str, required=True)
    return parser.parse_args()
    
    
if __name__ == "__main__":
    args = parse_options()
    junctions = create_junctions(43)
    rec = create_split_reconstruction(args.input_filename, junctions)
    stats = compute_halfsites_stats(rec)
    save_to_file(args.output_filename, stats)
    
