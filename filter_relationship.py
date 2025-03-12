# /usr/bin/python3.9
# filter_relationship.py
## read_annotation_and_relation adapted from https://github.com/katnastou/STRING-process-RE-output-multilabel/blob/78abc1ee6678e9f34d680f7a30365e1477fbf4c7/generate_output_for_RE_scoring_keep_brat_ids.py#L13
# input
## .outtsv.gz RE output file (or trigger input file) 
#pmid   e1      e2      Catalysis_of_ADP-ribosylation>  Catalysis_of_ADP-ribosylation<  Catalysis_of_SUMOylation>        Catalysis_of_SUMOylation<       Catalysis_of_acetylation>       Catalysis_of_acetylation<        Catalysis_of_acylation> Catalysis_of_acylation< Catalysis_of_deSUMOylation>      Catalysis_of_deSUMOylation<     Catalysis_of_deacetylation>     

## .ann file
# #911	AnnotatorNotes T264 T264|19362000|70|1|35338|35344|9031.ENSGALP00000016266
# #912	AnnotatorNotes T264 T264|19362000|70|1|35338|35344|9606.ENSP00000298159
# R1	Regulation Arg1:T16 Arg2:T15
# R2	Complex_formation Arg1:T26 Arg2:T25

# output
## same-org-outputs-with-eids.tsv
# 391	T8	T11	9606	ENSP00000303830	2	7	986	1001	9606	ENSP00000380432	2	7	1099	1105	Complex_formation	0.9997691512107849
# 391	T8	T9	9606	ENSP00000303830	2	7	986	1001	9606	ENSP00000380432	2	7	1016	1022	Complex_formation	0.752830445766449


import os
import sys
import shutil
import argparse
import gzip
import tarfile
import numpy as np

def read_annotation_and_relation(document_id, input_tar_gz_handle):
    doc_dict={}
    document_id = str(document_id)
    rooted_document_id = "./" + document_id

    document_id = int(document_id)
    doc_dict = {}
    rel_dict = {}
    
    tar_member_ann = input_tar_gz_handle.getmember(rooted_document_id + ".ann")
    f = input_tar_gz_handle.extractfile(tar_member_ann)
    lines = f.readlines()
    for line in lines:
        line = str(line.decode())
        if line[0] == "#":
            _ , e_id , e_notes = line.strip().split(" ")
            e_id, pmid, par_num, sent_num, start, end, string_id = e_notes.split("|")
            e_id, pmid, par_num, sent_num, start, end, string_id = str(e_id), int(pmid), int(par_num), int(sent_num), int(start), int(end), str(string_id)
            if pmid == document_id:
                if e_id not in doc_dict:
                    doc_dict[e_id] = {}
                    doc_dict[e_id]['par_num'] = par_num
                    doc_dict[e_id]['sent_num'] = sent_num
                    doc_dict[e_id]['start'] = start
                    doc_dict[e_id]['end'] = end
                    doc_dict[e_id]['string_id'] = [string_id]
                else:
                    doc_dict[e_id]['string_id'].append(string_id)
        elif line[0] == 'R':
            r_id, rest = line.strip().split('\t')
            rel_type, arg1, arg2 = rest.split()
            e1_id = int(arg1.split('T')[1])
            e2_id = int(arg2.split('T')[1])
            if e1_id < e2_id:
                k = (f'T{e1_id}', f'T{e2_id}')
                if k in rel_dict:
                    rel_dict[k].append((r_id, rel_type+'>'))
                else:
                    rel_dict[k] = [(r_id, rel_type+'>')]
            else:
                k = (f'T{e2_id}', f'T{e1_id}')
                if k in rel_dict:
                    rel_dict[k].append((r_id, rel_type+'<'))
                else:
                    rel_dict[k] = [(r_id, rel_type+'<')]
    return doc_dict, rel_dict


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tar_gz_ann_path" , type=str)
    parser.add_argument("output_gz_tsv_path" , type=str)
    parser.add_argument("filtered_relation_with_eid_filepath" , type=str)
    args = parser.parse_args()
    
    input_tar_gz_handle = tarfile.open(args.input_tar_gz_ann_path , "r:gz")
    pmid_doc_dict = {}
    pmid_rel_dict = {}
    with gzip.open(args.output_gz_tsv_path , "rt") as f:
        next(f)
        for l in f:
            document_id = l.split('\t')[0]
            doc_dict, rel_dict = read_annotation_and_relation(document_id, input_tar_gz_handle)
            pmid_doc_dict[document_id] = doc_dict
            pmid_rel_dict[document_id] = rel_dict
            
    with gzip.open(args.output_gz_tsv_path , "rt") as f:
        _ , _ , _ , *labels  = f.readline().strip().split("\t")
        for line in f:
            document_id , e1_id , e2_id , *probabilities_ = line.strip().split("\t")
            # probabilities=np.array(probabilities_).astype(float)
            # q = (1 - probabilities)
            # e1_id, e2_id = str(e1_id), str(e2_id)
            document_id = int(document_id)
            
            if (e1_id, e2_id) not in pmid_rel_dict[document_id]:
                continue
            
                
            for i in doc_dict[e1_id]['string_id']:
                for j in doc_dict[e2_id]['string_id']:
                    
                    taxid_e1, stringid_e1 = str(i).split(".",1) #stringid might have several .
                    taxid_e2, stringid_e2 = str(j).split(".",1)
                    
                    #print line if the organism is the same
                    if (taxid_e1 == taxid_e2):
                        #print line if the protein is not the same
                        if (str(i) != str(j)):
                            #print line if things are in the same paragraph only
                            if (pmid_doc_dict[document_id][e1_id]['par_num'] == pmid_doc_dict[document_id][e2_id]['par_num']):
                               
                                # write this line
                                # 391	T8	T11	9606	ENSP00000303830	2	7	986	1001	9606	ENSP00000380432	2	7	1099	1105	Complex_formation	0.9997691512107849
                                
                                ## can be multiple relations per pair
                                for (rel_id, rel_type) in pmid_rel_dict[document_id]: 
                                    label_idx = label.index(pmid_rel_dict[document_id][rel_type])
                                    relationship_score = 1-float(probabilities_[label_idx])
                                    
                                    with open(args.filtered_relation_with_eid_filepath, 'a') as wf:
                                        if rel_type[-1] == '<':
                                            e1_p = e2_id
                                            e2_p = e1_id
                                            tax1_p = taxid_e2
                                            str1_p = stringid_e2
                                            tax2_p = taxid_e1
                                            str2_p = stringid_e1
                                            rel_type_print = rel_type[:-1]
                                        elif rel_type[-1] == '>':
                                            rel_type_print = rel_type[:-1]
                                        else:
                                            rel_type_print = rel_type
                                            
                                        print(document_id, e1_p, e2_p, 
                                            tax1_p, str1_p, doc_dict[e1_p]['par_num'], doc_dict[e1_p]['sent_num'], doc_dict[e1_p]['start'], doc_dict[e1_p]['end'],
                                            tax2_p, str2_p, doc_dict[e2_p]['par_num'], doc_dict[e2_p]['sent_num'], doc_dict[e2_p]['start'], doc_dict[e2_p]['end'],
                                            rel_id, rel_type_print
                                            sep='\t', file=wf)