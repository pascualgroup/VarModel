#ifndef __malariamodel__DatabaseTypes__
#define __malariamodel__DatabaseTypes__

#include "zppdb.hpp"

using namespace zppdb;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/**
	\brief Row type for meta table.
*/
ZPPDB_DEFINE_ROW_TYPE(
	MetaRow,
	/**
		\brief meta-information key (`"parameters"`)
	*/
	((Text)(key))
	
	/**
		\brief value (parameter values in JSON format)
	*/
	((Text)(value))
)

/**
	\brief Type defining rows in `genes` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	GeneRow,
	
	/**
		\brief Unique identifier
	*/
	( (Integer)(geneId) )
	
	/**
		\brief Transmissibility parameter
	*/
	( (Real)(transmissibility) )
	
	/**
		\brief Rate of immunity loss to this gene
	*/
	( (Real)(immunityLossRate) )
	
	/**
		\brief source of genes, original pool = 0, recombination = 1, mutation = 2, immigration = 3
	*/
	( (Integer)(source) )
    /**
    \brief Functional gene or not
    */
    ( (Integer)(functionality) )
                      
    
)

/**
 \brief Type defining rows in `loci` table
 */
ZPPDB_DEFINE_ROW_TYPE(
    LociRow,
                      
   /**
   \brief Unique identifier of the gene
   */
   ( (Integer)(geneId) )

    /**
    Expression index of a locus identified by this row
    */
    ( (Integer)(alleleIndex) )
    
    /**
    Unique identifier of the allele at this index
    */
    ( (Integer)(alleleId) )
)

/**
	\brief Type defining rows in `strains` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	StrainRow,
	
	/**
		Unique identifier for strain
	*/
	( (Integer)(strainId) )

	/**
		Expression index of gene identified by this row
	*/
	( (Integer)(geneIndex) )

	/**
		Unique identifier of gene at this index
	*/
	( (Integer)(geneId) )
)

/**
	\brief Type defining rows in `hosts` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	HostRow,
	
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
 \brief Type defining rows in `followedHosts` table
*/
ZPPDB_DEFINE_ROW_TYPE(
     followedHostsRow,
     /**
      Sampling time
     */
     ( (Real)(time) )
     /**
     Unique identifier for host
     */
     ( (Integer)(hostId) )
                      
     /**
     Population id for the host
     */
     ( (Integer)(popId) )

     /**
     infectionId
     */
     ((Integer)(infectionId))
    
     /**
      strainId
      */
     ((Integer)(strainId))
      /**
       duration of infection
       */
     ((Real)(duration))
                      
)

/**
	\brief Type defining rows in `sampledHosts` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	SampledHostRow,
	
	/**
		Sampling time
	*/
	( (Real)(time) )
	
	/**
		Host ID
	*/
	( (Integer)(hostId) )
)

/**
	\brief Type defining rows in `sampledTransmissions` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	TransmissionRow,
	
	/**
		\brief ID of sampled transmission
	*/
	( (Integer)(transmissionId) )

	/**
		\brief ID of source host
	*/
	( (Integer)(sourceHostId) )

	/**
		\brief ID of target host
	*/
	( (Integer)(targetHostId) )
)

/**
	\brief Type defining rows in `sampledTransmissionStrains` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	TransmissionStrainRow,
	/**
        \brief time
    */
     ((Real)(time))
	/**
		\brief ID of sampled transmission
	*/
	( (Integer)(transmissionId) )

	/**
		\brief ID of transmitted strain
	*/
	( (Integer)(strainId) )
)

/**
	\brief Type defining rows in `sampledTransmissionInfections` table
*/
ZPPDB_DEFINE_ROW_TYPE(
	TransmissionInfectionRow,
	
	/**
		\brief Time of infection sampling
	*/
	( (Integer)(transmissionId) )
	
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
		\brief Index of gene within strain
	*/
	( (Integer)(geneIndex) )
	
	/**
		Whether or not the gene is currently active
	*/
	( (Integer)(active) )
)

/**
	\brief Type defining rows in `sampledTransmissionImmunity and
	`sampledTransmissionClinicalImmunity` tables.
*/
ZPPDB_DEFINE_ROW_TYPE(
	TransmissionImmunityRow,
	
	/**
		Transmission ID
	*/
	( (Integer)(transmissionId) )
	
	/**
		Host ID
	*/
	( (Integer)(hostId) )
	
	/**
		Gene ID
	*/
	( (Integer)(geneId) )
	
	/**
		Loss rate of immunity
	*/
	( (Real)(lossRate) )
)

/**
	\brief Type defining rows in `sampledHostInfections` and
	`sampledTransmissionInfections` tables
*/
ZPPDB_DEFINE_ROW_TYPE(
	InfectionRow,
	
	/**
		\brief Time of infection sampling
	*/
	( (Real)(time) )
	
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
		\brief Index of gene within strain
	*/
	( (Integer)(geneIndex) )
	
	/**
		Whether or not the gene is currently active
	*/
	( (Integer)(active) )

   /**
   \brief the microsat id of the current infection
   */
   ((Integer) (msID))
)

/**
	\brief Type defining rows in `sampledHostImmunity`,
	`sampledHostClinicalImmunity`, `sampledTransmissionImmunity, and
	`sampledTransmissionClinicalImmunity` columns.
*/
ZPPDB_DEFINE_ROW_TYPE(
    InfectionDurationRow,
    /**
     time of infection sampled started
     */
     ((Real)(time))
     /**
      duration of infection
      */
     ((Real)(duration))
    /**
     hostId
     */
     ((Integer)(hostId))
      /**
      infectionId
      */
     ((Integer)(infectionId))
)

ZPPDB_DEFINE_ROW_TYPE(
	ImmunityRow,
	
	/**
		Time of immunity sampling
	*/
	( (Real)(time) )
	
	/**
		Host ID
	*/
	( (Integer)(hostId) )
	
	/**
		Gene ID
	*/
	( (Integer)(geneId) )
	
	/**
		Loss rate of immunity
	*/
	( (Real)(lossRate) )
)

/**
 \brief Type defining rows `alleleImmunity` columns.
 */
ZPPDB_DEFINE_ROW_TYPE(
    AlleleImmunityRow,
    /** Time of immunity occurred
     */
    ((Real)(time))
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
)

/**
 \brief Type defining rows `recordEIR` columns.
 */
ZPPDB_DEFINE_ROW_TYPE(
   recordEIRRow,
   /** Time of bit event occurred
   */
   ((Real)(time))
   /**
    Whether the bit is infectious
   */
   ((Integer)(infectious))
)

#endif // #define DOXYGEN_SHOULD_SKIP_THIS

#endif // #ifndef __malariamodel__DatabaseTypes__
