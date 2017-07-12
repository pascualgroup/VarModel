# Checkpointing notes

## Saving to checkpoint

Set these keys in the parameters file to indicate where and how often (in simulation time) checkpoints should be saved:

```json
{	
    "checkpointSaveFilename" : "checkpoint.sqlite",
    "checkpointPeriod" : 360
}
```

The `Simulation::saveCheckpoint()` function is called every `checkpointPeriod` simulation time units; it does the following things in its call tree:

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
    `infections` table: strain, microsat, expression index, etc. for each infection
    `expressionOrder` table: expression order for each active infection
    `alleleImmunity` table: allele-based immunity (if `useAlleleImmunity == true`)
    `immunity` table: gene-based immunity (if `useAlleleImmunity == false`)


## Loading from checkpoint

Use this key to indicate what checkpoint to load from:

```json
{
	"checkpointLoadFilename" : "checkpoint.sqlite"
}
```

The `Simulation::loadCheckpoint()` function loads the simulation from the saved checkpoint by doing the folloiwng:

* Load genes and microsats: these functions get called with regular gene and microsat data:
    * `loadAlleleCounts()`
    * `loadGenes()`
    * `loadStrains()`
    * `loadPopulations()`
        * Still To Finish:: hosts, infections, immune history
