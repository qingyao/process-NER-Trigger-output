# process-NER-Trigger-output

Re-write scripts to process the NER-based trigger detection output to merge with its input to generate evidence viewer.

## filter_relationship.py: keep interactions from same organism. keep segment for documents with interactions. Filter the original regulatome prediction based on the extracted relations (in .ann files)

---
Input

1. .outtsv.gz RE output file (or trigger input file) 
*Here restriction e1_id < e2_id*
```#pmid   e1      e2      Catalysis_of_ADP-ribosylation>  Catalysis_of_ADP-ribosylation<  Catalysis_of_SUMOylation>        Catalysis_of_SUMOylation<       Catalysis_of_acetylation>       Catalysis_of_acetylation<        Catalysis_of_acylation> Catalysis_of_acylation< Catalysis_of_deSUMOylation>      Catalysis_of_deSUMOylation<     Catalysis_of_deacetylation>     Catalysis_of_deacetylation<      Catalysis_of_deacylation>       Catalysis_of_deacylation<       Catalysis_of_deglycosylation>    Catalysis_of_deglycosylation<   Catalysis_of_demethylation>     Catalysis_of_demethylation<      Catalysis_of_deneddylation>     Catalysis_of_deneddylation<     Catalysis_of_depalmitoylation>   Catalysis_of_depalmitoylation<  Catalysis_of_dephosphorylation> Catalysis_of_dephosphorylation<  Catalysis_of_deubiquitination>  Catalysis_of_deubiquitination<  Catalysis_of_farnesylation>      Catalysis_of_farnesylation<     Catalysis_of_geranylgeranylation>       Catalysis_of_geranylgeranylation<        Catalysis_of_glycosylation>     Catalysis_of_glycosylation<     Catalysis_of_lipidation> Catalysis_of_lipidation<        Catalysis_of_methylation>       Catalysis_of_methylation<        Catalysis_of_neddylation>       Catalysis_of_neddylation<       Catalysis_of_other_small_molecule_conjugation_or_removal>        Catalysis_of_other_small_molecule_conjugation_or_removal<        Catalysis_of_palmitoylation>    Catalysis_of_palmitoylation<    Catalysis_of_phosphoryl_group_conjugation_or_removal>    Catalysis_of_phosphoryl_group_conjugation_or_removal<   Catalysis_of_phosphorylation>    Catalysis_of_phosphorylation<   Catalysis_of_posttranslational_modification>     Catalysis_of_posttranslational_modification<    Catalysis_of_prenylation>       Catalysis_of_prenylation<        Catalysis_of_small_molecule_removal>    Catalysis_of_small_molecule_removal<     Catalysis_of_small_protein_conjugation> Catalysis_of_small_protein_conjugation< Catalysis_of_small_protein_conjugation_or_removal>       Catalysis_of_small_protein_conjugation_or_removal<       Catalysis_of_small_protein_removal>     Catalysis_of_small_protein_removal<     Catalysis_of_ubiquitination>     Catalysis_of_ubiquitination<    Complex_formation       Negative_regulation>     Negative_regulation<    Other_catalysis_of_small_molecule_conjugation>  Other_catalysis_of_small_molecule_conjugation<   Other_catalysis_of_small_molecule_removal>      Other_catalysis_of_small_molecule_removal<       Other_catalysis_of_small_protein_conjugation>   Other_catalysis_of_small_protein_conjugation<    Other_catalysis_of_small_protein_removal>       Other_catalysis_of_small_protein_removal<        Positive_regulation>    Positive_regulation<    Regulation>      Regulation<     Regulation_of_degradation>      Regulation_of_degradation<      Regulation_of_gene_expression>   Regulation_of_gene_expression<  Regulation_of_transcription>    Regulation_of_transcription<     Regulation_of_translation>      Regulation_of_translation<
10220000        T1      T2      6.883724146256043e-10   2.239807672665961e-09   6.661561080534284e-09    2.324611170223534e-09   8.242782456591158e-09   1.4324108565944016e-09  2.9225960540557594e-10   1.8656432301811243e-11  1.9222831659782003e-10  1.7725195768612811e-13  9.054130223340451e-10    1.223164625141493e-10   2.3418893560500642e-11  1.7446131616184052e-13  1.1755466391982772e-12   1.7187781146260628e-13  2.665250908862049e-09   1.052563591841249e-09   1.3195912140773203e-09   2.1970408131677388e-11  2.699481749246502e-10   1.9353476341053139e-13  7.049215433596601e-09    3.369109435880091e-09   2.9430375914074602e-09  1.2171097463209435e-09  4.5792362957097765e-11   1.2197549827813736e-12  7.607404151066532e-10   2.457000054800762e-11   2.5099109457471513e-09   1.9317169391852573e-11  1.231971664440723e-12   4.819166333352998e-12   6.492951953873671e-09    5.9007922992293516e-09  2.040846824868936e-09   1.8617702865486585e-13  7.269299051593237e-11    1.3448735126844052e-13  2.476674643148158e-09   1.935941129627139e-10   5.1321976182738815e-12   6.88251261538042e-11    4.768575223579319e-08   2.955797917536529e-08   6.8638912331664415e-09   5.177580142401439e-09   2.6234370578692712e-11  1.5313218629329356e-13  4.726685189082591e-12    1.3973993132337242e-13  2.795989217929673e-09   1.4143436422031641e-09  6.17160706220532e-12     1.3554584971792183e-12  1.2457616318725662e-12  1.7720125768203726e-13  6.759631787645048e-08    2.376028795936236e-08   2.4204866804211633e-06  6.550284865625144e-07   3.305437985545723e-07    3.836314821370479e-09   6.906191764299852e-12   1.3211707118945815e-12  5.7963061433907725e-12   5.82410300542513e-11    2.663014665553465e-13   1.1778515445967641e-12  9.15179715864356e-13     9.248662991012679e-07   2.7245184242019604e-07  8.081080977717647e-07   1.798310194089936e-07    8.439047149977341e-08   2.0198971384388642e-08  2.2066683413868304e-07  3.6030332495329276e-08   2.441299500333116e-07   8.297379139321492e-08   3.6415401249989543e-10  5.238624204567666e-12
```

2. .ann file
* All AnnotatorNotes of same entity are used, as long as the pair are from same organism. 
* All relations are taken to find the column in .outtsv.gz
*Here in relation: Arg1 can have larger id than Arg2, but relation has no direction.*
```
T263	Protein 34660 34667	cofilin
#906	AnnotatorNotes T263 T263|19362000|68|1|34660|34666|284812.P78929
#907	AnnotatorNotes T263 T263|19362000|68|1|34660|34666|9031.ENSGALP00000016266
#908	AnnotatorNotes T263 T263|19362000|68|1|34660|34666|9606.ENSP00000298159
#909	AnnotatorNotes T263 T263|19362000|68|1|34660|34666|9606.ENSP00000432660
T264	Protein 35338 35345	cofilin
#910	AnnotatorNotes T264 T264|19362000|70|1|35338|35344|284812.P78929
#911	AnnotatorNotes T264 T264|19362000|70|1|35338|35344|9031.ENSGALP00000016266
#912	AnnotatorNotes T264 T264|19362000|70|1|35338|35344|9606.ENSP00000298159
R1	Regulation Arg1:T16 Arg2:T15
R2	Complex_formation Arg1:T26 Arg2:T25
R3	Complex_formation Arg1:T116 Arg2:T115
R4	Regulation Arg1:T119 Arg2:T118
```

---
Output
1. same-org-outputs-with-eids.tsv
*Here also e2_id>e1_id because*
`awk  '{match($2, /T([0-9]+)/, a);match($3, /T([0-9]+)/, b);if (a[1]>b[1]) {print}}' results/same-org-outputs-with-eids.tsv|head` *returns nothing and also for Complex_formation, it didn't matter*
```
*but now revert this as final file can have e2_id < e1_id.*
391	T8	T11	9606	ENSP00000303830	2	7	986	1001	9606	ENSP00000380432	2	7	1099	1105	Complex_formation	0.9997691512107849
391	T8	T9	9606	ENSP00000303830	2	7	986	1001	9606	ENSP00000380432	2	7	1016	1022	Complex_formation	0.752830445766449
401	T1	T2	9606	ENSP00000317159	1	1	21	33	9606	ENSP00000307786	1	1	39	50	Complex_formation	0.9998788833618164
```

2. segments-in-docs-with-interactions-with-eids.tsv
Maybe not needed, for now. 
```
391	1	1	0	47
391	2	1	108	254
391	1	2	48	106
391	2	2	255	420
```


## merge_pair_trigger_to_evidence.py: merge input (.ann, .txt file) and output (trigger detected) from trigger pipeline into evidence viewer format.

---
Input
1. .ann file 
```
T263	Protein 34660 34667	cofilin
#906	AnnotatorNotes T263 T263|19362000|68|1|34660|34666|284812.P78929
#907	AnnotatorNotes T263 T263|19362000|68|1|34660|34666|9031.ENSGALP00000016266
#908	AnnotatorNotes T263 T263|19362000|68|1|34660|34666|9606.ENSP00000298159
#909	AnnotatorNotes T263 T263|19362000|68|1|34660|34666|9606.ENSP00000432660
T264	Protein 35338 35345	cofilin
#910	AnnotatorNotes T264 T264|19362000|70|1|35338|35344|284812.P78929
#911	AnnotatorNotes T264 T264|19362000|70|1|35338|35344|9031.ENSGALP00000016266
#912	AnnotatorNotes T264 T264|19362000|70|1|35338|35344|9606.ENSP00000298159
R1	Regulation Arg1:T16 Arg2:T15
R2	Complex_formation Arg1:T26 Arg2:T25
R3	Complex_formation Arg1:T116 Arg2:T115
R4	Regulation Arg1:T119 Arg2:T118
```

2. same-org-outputs-with-eids.tsv from previous step `filter_relationship.py`
```
391	T8	T11	9606	ENSP00000303830	2	7	986	1001	9606	ENSP00000380432	2	7	1099	1105	Complex_formation	0.9997691512107849
391	T8	T9	9606	ENSP00000303830	2	7	986	1001	9606	ENSP00000380432	2	7	1016	1022	Complex_formation	0.752830445766449
401	T1	T2	9606	ENSP00000317159	1	1	21	33	9606	ENSP00000307786	1	1	39	50	Complex_formation	0.9998788833618164
```
...

3. .outtsv.gz trigger_output file 
```
19362000,R5,T131,T130,Complex_formation,8.473395,15730,15735,binds
19362000,R6,T132,T130,Complex_formation,8.313196,15730,15735,binds
19362000,R7,T145,T144,Complex_formation,4.7857723,17838,17846,subunits
```


---
Output

1. .tsv evidence viewer file (selected-pairs-coordinates-with-probabilities-same-org-AB-BA-with-triggers.tsv)
* it used to have segment_start,segment_end, but for now not needed.
* instead add relation_type
*Here entity1 can appear later in text than entity2, just so that the relation will be clean from directionality, e.g. instead of Positive_regulation> and Positive_regulation<, having only Positive_regulation*
#PMID  entity1_start     entity1_end entity1_taxid      entity1_id  entity2_start     entity2_end entity2_taxid     entity2_id   relationship_type   relationship_score  triggerstart    triggerend  triggermatch    triggerscore
```

35500116	-(31256)-	-(31460)-	31396	31399	48	AA314_04725	31447	31450	48	AA314_09029	0.9998824596405029	31452	31458  +(Complex_formation)+	 binding	0.673339307584731
```

