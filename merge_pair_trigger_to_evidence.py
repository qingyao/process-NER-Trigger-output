import csv
import argparse
import os, sys
import glob, gzip
from collections import Counter, defaultdict

# input1 (** match key, * info column to extract)
# **pmid, e1_id, e2_id, *tax1_id, *str1_id, par1_num, sent1_num, *e1_start, *e1_end, *tax2_id, *str2_id, par2_num, sent2_num, *e2_start, *e2_end, **rel_id, rel_type, *rel_score
# 391	T8	T11	9606	ENSP00000303830	2	7	986	1001	9606	ENSP00000380432	2	7	1099	1105	+(R5)+   Complex_formation	0.9997691512107849

# input2 (** match key, * info column to extract)
# **pmid,**rel_id,*e1_id,*e2_id,*relation_type,*confidence_score,*trigger_bgn_offset,*trigger_end_offset,*trigger_text
# 19362000,R5,T131,T130,Complex_formation,8.473395,15730,15735,binds
# 19362000,R6,T132,T130,Complex_formation,8.313196,15730,15735,binds
# 19362000,R7,T145,T144,Complex_formation,4.7857723,17838,17846,subunits

# output1
#PMID  entity1_start     entity1_end entity1_taxid      entity1_id  entity2_start     entity2_end entity2_taxid     entity2_id   relationship_type   relationship_score  triggerstart    triggerend  triggermatch    trigger_id   triggerscore

# output2 stats
# for each relation, how many has trigger words at all, how often multiple trigger words, how the trigger scores differ.


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_same_org_eids_path" , type=str)
    parser.add_argument("--in_trigger_detection_tsvgz_path" , type=str)
    parser.add_argument("--out_evidence_viewer_path" , type=str)
    parser.add_argument("--out_log_path" , type=str)
    parser.add_argument("--out_stats_path" , type=str)
    args = parser.parse_args()
    
    if not os.path.isfile(args.out_evidence_viewer_path):
        with gzip.open(args.out_evidence_viewer_path, 'wt') as wf:
            print('pmid', 'entity1_start', 'entity1_end', 'entity1_taxid', 'entity1_id', 'entity2_start', 'entity2_end', 'entity2_taxid', 'entity2_id', 'relationship_type', 'relationship_score', 'triggerstart', 'triggerend', 'triggermatch', 'trigger_id', 'triggerscore,', sep= '\t', file = wf)
        
    if os.path.isfile(args.out_log_path):
        os.remove(args.out_log_path)
        
    eid_lookup_dict = {}
    with gzip.open(args.in_same_org_eids_path, 'rt', newline='', encoding='utf-8') as gz_f:
        reader = csv.DictReader(gz_f, delimiter = '\t')
        for row in reader:
            pmid = row['pmid']
            rel_id = row['rel_id']
            tax1_id = row['tax1_id']
            str1_id = row['str1_id']
            e1_start = row['e1_start']
            e1_end = row['e1_end']
            tax2_id = row['tax2_id']
            str2_id = row['str2_id']
            e2_start = row['e2_start']
            e2_end = row['e2_end']
            rel_score = row['rel_score']
            eid_lookup_dict[pmid+rel_id] = [tax1_id, str1_id, e1_start, e1_end, tax2_id, str2_id, e2_start, e2_end, rel_score]
    
    trigger_count = {}
    trigger_scores = defaultdict(list)
    rel_type_to_pmid_relid = defaultdict(list)
    for filepath in glob.glob(args.in_trigger_detection_tsvgz_path):
        with gzip.open(filepath, 'rt', newline='', encoding='utf-8') as gz_f:
            
            reader = csv.DictReader(gz_f)
            
            for row in reader:
                pmid = row['pmid']
                rel_id = row['rel_id']
                e1_id = row['e1_id']
                e2_id = row['e2_id']
                rel_type = row['relation_type']
                trig_st = int(row['trigger_bgn_offset'])
                trig_end = int(row['trigger_end_offset'])
                trigger_score = float(row['confidence_score'])
                trigger_text = row['trigger_text']
                
                ## may have been excluded if e1_id and e2_id is not in the same paragraph
                ## only continue if the pair exists
                if pmid+rel_id in eid_lookup_dict: 
                    
                    # handle new line and change offsets
                    if trigger_text.count('\n') == 1:
                        rm_pos = trigger_text.index('\n')
                        if rm_pos == 0:
                            trig_st += 1
                        elif rm_pos == len(trigger_text) -1:
                            trig_end -= 1
                        else:
                            ## log
                            with open(args.out_log_path, 'a') as wf:
                                print(pmid, rel_id, 'newline in middle of trigger', file = wf)
                            continue
                    elif trigger_text.count('\n') > 1:
                        ## multiple new lines
                        ## log
                        with open(args.out_log_path, 'a') as wf:
                            print(pmid, rel_id, 'trigger contains multiple newlines', file = wf)
                        continue
                    
                    if trig_st == trig_end:
                        ## log
                        with open(args.out_log_path, 'a') as wf:
                            print(pmid, rel_id, 'trigger contains only newline', file = wf)
                        continue
                    
                    rel_type_to_pmid_relid[rel_type].append(pmid+rel_id)
                    trigger_scores[pmid+rel_id].append(trigger_score) ## reduce score precision to save stats write-out space
                    
                    if pmid+rel_id in trigger_count:
                        trigger_count[pmid+rel_id] += 1
                    else:
                        trigger_count[pmid+rel_id] = 1
                
                    with gzip.open(args.out_evidence_viewer_path, 'at') as wf:
                        tax1_id, str1_id, e1_start, e1_end, tax2_id, str2_id, e2_start, e2_end, rel_score = eid_lookup_dict[pmid+rel_id]
                        
                        print(pmid, e1_start, e1_end, tax1_id, str1_id, e2_start, e2_end, tax2_id, str2_id, rel_type, rel_score, trig_st, trig_end, trigger_text, f'Tr_{trigger_count[pmid+rel_id]}', trigger_score, sep = '\t', file = wf)
                    
    ## stats
    
    with open(args.out_stats_path, 'w') as wf:
        
        for rel_type, rel_info in rel_type_to_pmid_relid.items():
            print(rel_type, file = wf)
            print('Trigger scores when multiple triggers are found:', file = wf)
            tw_counts = []
            for pmid_rel_id in rel_info:
                if pmid_rel_id in trigger_count:
                    tw_counts.append(trigger_count[pmid_rel_id])
                    if trigger_count[pmid_rel_id] > 1:
                        print(','.join([f'{i:.2f}' for i in sorted(trigger_scores[pmid_rel_id])]), file=wf)
                        
                else:
                    tw_counts.append(0)
            
            print('Trigger counts per relationship:', file = wf)
            for tw_num, counts in Counter(tw_counts).items():
                print(tw_num, counts, sep = ':', file = wf)

            print('-----------------------------', file = wf)
                
    