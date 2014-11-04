# Introduction

Analysis of Structural Variation events has been done in the 1000 Genomes project and similar project as GoNL on population scale. These projects provided good insights into the characteristics of deletion and insertion events, aswell a little collection of translocations.

The toolset we provide in `pySVTools` is based on the theory and experiences gained from both projects. The implementation is solely based on conversations with participants in the projects. Ideas were adapted and prototypes were accessed on real tumor datasets. 

By releasing `pySVTools` we aim at refinement of current tools and possibly helping others.

# Concept

## Define overlap

Conventional methods use the reciprocal overlap to determine whether a genomic event is the same or not.
This approach works most of deletion events, by given any arbitrary number of bases overlap. This approach however doesn't identify 
whether an overlapping event is describing the same biological consequence. The `overlap` is solely based on position-wise matching.

## Intersection calls

The power of intersection is the ability to find shared denominators. In Structural Variation analysis this is a common task if one is working in the field of tumor genetics. We have worked out the following concept:

 * A SVcall should have `overlapping` centerpoints, default setting with `-f 100`, denoting we flank the centerpoint by 100bp. In technical terms, we do a reciprocal overlap on the flanked centerpoint definition. In which should produce a very sensitive merge operation reducing false positives.
 * A SVcall should have similar size to be called catagorized into the same catagory.
 * A SVcall appears in at least 2 target samples (tumor sequences)
 * Control samples can be supplied for filtering purposes. (keeping the all records in, just an annotation now)

# Interpretation of results

A sample report of the intersected results:

|ChrA|	|ChrApos|	|ChrB|	|ChrBpos|	|SVTYPE|	|DP|	|Size|	|FA.commonhits.vcf|	|size|	|530.commonhits.vcf|	|size|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|chr10|	|2565649|	|chr10|	|2565806|	|INS|	|38|	|157|	|chr10:2565762-2565886|	|124|	|chr10:2565649-2565806|	|157|
|chr10|	|3095759|	|chr10|	|3095930|	|DEL|	|26|	|171|	|chr10:3095199-3096437|	|1238|	|chr10:3095129-3096551|	|1422|
|chr10|	|3690412|	|chr10|	|3691814|	|INS|	|61|	|1402|	|chr10:3690412-3691814|	|1402|	|	|	|	|
|chr10|	|4290075|	|chr10|	|4291663|	|DEL|	|73|	|1588|	|chr10:4290103-4291692|	|1589|	|chr10:4290109-4291673|	|1564|

Columns explained:

 1. Chromosome of breakpoint 1
 1. Chromosomal position of breakpoint 1
 1. Chromosome of breakpoint 2
 1. Chromosomal position of breakpoint 2
 1. Annotated Structural Variation Type, from: DEL, INS, CTX, INV, ITX
 1. Estimated Size of the event

Sample specific columns

 1. Chromosomal position of both ends on this sample
 1. Estimated size of this event in this sample

## 
