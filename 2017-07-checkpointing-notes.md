# Checkpointing notes

## Saving to checkpoint

Set these keys in the parameters file to indicate where and how often (in simulation time) checkpoints should be saved:

```json
{	
    "checkpointSaveFilename" : "checkpoint.sqlite",
    "checkpointPeriod" : 360
}

```

Checkpoint save call tree:

* `writeMetaToCheckpoint()`: saves one row to `meta` table
    * `time`
    * `parameters` (as JSON)
    * `nextHostId`
    * `nextStrainId`
* `writeGenesAndStrainsToCheckpoint()`
    * `alleleCounts` table: number of alleles at each locus
    * `genes` table: gene metadata
    * `geneAlleles` table: allele composition of each gene
    * `strains` table: gene composition of each strain
* `writeMicrosatsToCheckpoint()`
    * `microsatAlleleCounts` table: number of alleles at each microsat locus
    * `microsats` table: metadata on microsats
    * `microsatAlleles` table: allele composition of microsats
* `writePopulationsToCheckpoint()`
    `hosts` table: metadata for hosts
    `infections` table: gene, expression index, etc. for each infection
    `expressionOrder` table: expression order for each active infection
    `alleleImmunity` table: allele-based immunity (if `useAlleleImmunity == true`)
    `immunity` table: gene-based immunity (if `useAlleleImmunity == false`)
