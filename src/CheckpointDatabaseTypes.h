#ifndef __malariamodel__CheckpointDatabaseTypes__
#define __malariamodel__CheckpointDatabaseTypes__

#include "zppdb.hpp"

using namespace zppdb;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/**
    \brief Row type for checkpoint meta table.
*/
ZPPDB_DEFINE_ROW_TYPE(
    CheckpointMetaRow,
    
    /**
        \brief Time checkpoint was written out.
    */
    ((Real)(time))
    
    /**
        \brief Parameter values in JSON format.
    */
    ((Text)(parameters))
    
    /**
        \brief Simulation::nextHostId.
    */
    ((Integer)(nextHostId))
    
    /**
        \brief Simulation::nextStrainId.
    */
    ((Integer)(nextStrainId))
)

/**
	\brief Type defining rows in checkpoint `hosts` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	CheckpointHostRow,
	
	/**
		Unique identifier for host
	*/
	( (Integer)(hostId) )

    /**
    Population id for the host
    */
    ( (Integer)(popId) )
	/**
		Host birth time
	*/
	( (Real)(birthTime) )

	/**
		Host death time
	*/
	( (Real)(deathTime) )
)

/**
	\brief Type defining rows in checkpoint `infections` table.
*/
ZPPDB_DEFINE_ROW_TYPE(
	CheckpointInfectionRow,
	
	/**
		\brief Host ID
	*/
	( (Integer)(hostId) )
	
	/**
		\brief Infection ID
	*/
	( (Integer)(infectionId) )
	
	/**
		\brief Strain ID
	*/
	( (Integer)(strainId) )

   /**
   \brief The microsatellite id of the current infection
   */
   ((Integer) (msId))
	
	/**
		\brief Index of gene within strain
	*/
	( (Integer)(geneIndex) )
	
	/**
		\brief Whether or not the gene is currently active
	*/
	( (Integer)(active) )
    
    /**
        \brief The time of the last transition
    */
    ( (Real)(transitionTime) )
    
    /**
        \brief The initial time of infection
    */
    ( (Real)(initialTime) )
)

/**
 \brief Type defining rows in the checkpoint `expressionOrder` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	CheckpointExpressionOrderRow,
	
	/**
		\brief Host ID
	*/
	( (Integer)(hostId) )
	
	/**
		\brief Infection ID
	*/
	( (Integer)(infectionId) )
	
	/**
		\brief Expression order
	*/
	( (Integer)(expressionOrder) )
	
	/**
		\brief Index of gene within strain
	*/
	( (Integer)(geneIndex) )
)

/**
 \brief Type defining rows in checkpoint `immunity` table.
*/
ZPPDB_DEFINE_ROW_TYPE(
	CheckpointImmunityRow,
	
	/**
		Host ID
	*/
	( (Integer)(hostId) )
	
	/**
		Gene ID
	*/
	( (Integer)(geneId) )
)

/**
 \brief Type defining rows in checkpoint `alleleImmunity` table.
*/
ZPPDB_DEFINE_ROW_TYPE(
    CheckpointAlleleImmunityRow,
    
    /**
    Host ID
    */
    ((Integer)(hostId))
    
    /**
     Index of the locus in a gene
     */
    ((Integer)(locusIndex))
    
    /**
     Allele id of the locus
     */
    ((Integer)(alleleId))
    
    /**
    Immunity level (number of exposures)
    */
    ((Integer)(level))
)

/**
 \brief Type defining rows in checkpoint `allele_counts` table (number of mutated alleles)
*/
ZPPDB_DEFINE_ROW_TYPE(
    CheckpointAlleleCountRow,
    
    /**
     Index of the locus in a gene
     */
    ((Integer)(locusIndex))
    
    /**
     Number of alleles at that locus
     */
    ((Integer)(alleleCount))
)

#endif // #define DOXYGEN_SHOULD_SKIP_THIS

#endif // #ifndef __malariamodel__CheckpointDatabaseTypes__
