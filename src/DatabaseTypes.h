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
		\brief Rate of loss of "clinical immunity"
	*/
	( (Real)(clinicalImmunityLossRate) )
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
		Host birth time
	*/
	( (Real)(birthTime) )

	/**
		Host death time
	*/
	( (Real)(deathTime) )
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
		\brief Time of sampled transmission
	*/
	( (Real)(time) )

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
)

/**
	\brief Type defining rows in `sampledHostImmunity`,
	`sampledHostClinicalImmunity`, `sampledTransmissionImmunity, and
	`sampledTransmissionClinicalImmunity` columns.
*/
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

#endif // #define DOXYGEN_SHOULD_SKIP_THIS

#endif // #ifndef __malariamodel__DatabaseTypes__
