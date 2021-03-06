#ifndef __malariamodel__SimParameters__
#define __malariamodel__SimParameters__

#include "zppjson.hpp"

using namespace zppjson;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/**
	\brief Type for defining sinusoidal variable (e.g., seasonal biting rate).
	
	The value of the variable is equal to
	`mean * (1.0 + relativeAmplitude * sin(2 * pi * ((t / period) - phase)))`
*/
ZPPJSON_DEFINE_TYPE(
	Sinusoid,
	
	/**
		\brief Mean value of seasonal variable.
	*/
	( (Double)(mean) )
	
	/**
		\brief Amplitude of seasonal variable, as a fraction of mean.
	*/
	( (Double)(relativeAmplitude) )
	
	/**
		\brief Period of seasonal variable, in simulation time units.
	*/
	( (Double)(period) )
	
	/**
		\brief Phase of the seasonal variable as a fraction, between 0 and 1,
		of the period.
	*/
	( (Double)(phase) )
)

double evaluateSinusoid(Sinusoid & s, double t);


/**
	\brief Type representing discretized distribution.
*/
ZPPJSON_DEFINE_TYPE(
	DiscretizedPDF,
	
	/**
		\brief Lower bound of the distribution's support.
	*/
	( (Double)(x0) )
	
	/**
		\brief Discretization of PDF, vector of intervals of dx.
	*/
	( (Array<Double>)(dx) )
	
	/**
		\brief Relative probability of each discrete segment of PDF.
		
		Probabilities map to `[x0 + i * dx, x0 + (i + 1) * dx]`.
	*/
	( (Array<Double>)(pdf) )
)

/**
	\brief Type defining parameters for a single population.
*/
ZPPJSON_DEFINE_TYPE(
	PopulationParameters,
	
	/**
		\brief Number of individuals in the population.
	*/
	( (Int64)(size) )
	
	/**
		\brief Number of hosts to sample from population.
	*/
	( (Int64)(sampleSize) )
	
    /**
        \brief The sample size specified above means sampling at least n number of MOI
         or n number of infected. If True, MOI 1 of sampleSize will be sampled; 
        otherwise, sampleSize of infected will be sampled;
     */
    ((Bool)(moi1))
	/**
		\brief Number of initial infections in population.
	*/
	( (Int64)(nInitialInfections) )

	/**
		\brief Sinusoidal "paired" biting rate per host.
		
		The probability per unit time that a host will be bitten in order to
		serve as the source of a transmission event.
	*/
	( (Sinusoid)(bitingRate) )

    /**
        \brief biting rate distribution, specified as different rate per month
    */
    ( (Array<Double>)(monthlyBitingRateDistribution) )

	/**
		\brief Immigration rate, in number of random infection events per unit
		time (per population, not per capita).
	*/
	( (Double)(immigrationRate) )
	
	/**
		\brief Probability that an immigration event will include new genes
		added to the pool.
	*/
	( (Double)(pImmigrationIncludesNewGenes) )
	
	/**
		\brief Number of new genes added to the pool during immigration event
		that includes new genes.
	*/
	( (Int64)(nImmigrationNewGenes) )
	
	/**
		\brief x-position of population in 2D space.
	*/
	( (Double)(x) )

	/**
		\brief y-position of population in 2D space.
	*/
	( (Double)(y) )
	
	/**
		\brief Effective "distance" of a population to itself
	*/
	( (Double)(selfDistance) )
)

/**
	\brief Type defining within-host parameters.
	
	These parameters are a starting point for within-host model development;
	the actual rules will likely be different from these, and these parameters
	will need to evolve accordingly.
*/
ZPPJSON_DEFINE_TYPE(
	WithinHostParameters,
	
	/**
		\brief Whether or not to divide transmissibility by #infections ("cost")
	found not used
	( (Bool)(coinfectionsReduceTransmission) )*/

    /**
    \brief whether immunity should be stored as allele level instead of gene level. if true, storeClinicalImmunity must be false
                     
    activationRate = activationRateConstant * nActiveInfections^activationRatePower
                     */
    ( (Bool)(useAlleleImmunity) )

	/**
		\brief The constant controlling the relationship between # active infections
		and activation rate.
		
		activationRate = activationRateConstant * nActiveInfections^activationRatePower
	*/
	( (Double)(activationRateConstant) )

	/**
		\brief The power controlling the relationship between # active infections
		and activation rate.
	*/
	( (Double)(activationRatePower) )

	/**
		\brief The constant controlling the relationship between # active infections
		and activation rate.
		
		deactivationRate = deactivationRateConstant * nActiveInfections^deactivationRatePower
	*/
	( (Double)(deactivationRateConstant) )

	/**
		\brief The power controlling the relationship between # active infections
		and deactivation rate.
	*/
	( (Double)(deactivationRatePower) )

	/**
		\brief The constant controlling the relationship between # active infections
		and clearance rate when not immune.
		
		If inactive, clearanceRate = 0.
		If active and not immune:
		clearanceRate = clearanceRateConstantNotImmune * nActiveInfections^clearanceRatePower
	*/
	( (Double)(clearanceRateConstantNotImmune) )

	/**
		\brief The power controlling the relationship between # active infections
		and clearance rate when immune.
		
		If inactive, clearanceRate = 0.
		If active and immune:
		clearanceRate = clearanceRateConstantImmune * nActiveInfections^clearanceRatePower
	*/
	( (Double)(clearanceRateConstantImmune) )

	/**
		\brief The power controlling the relationship between # active infections
		and clearance rate.
	*/
	( (Double)(clearanceRatePower) )
                    
    /**
        \brief The number of times required to gain general immunity
     */
    ((Double)(infectionTimesToImmune))
    /**
    \brief general immunity infection duration line fitting with selection, a, b, c, d
    */
    ( (Array<Double>)(generalImmunityParams) )
    
    /**
     \brief maximum MOI a host can get
     */
    ((Int64)(maxMOI))
)

/**
	\brief Type defining parameters governing genes.
*/
ZPPJSON_DEFINE_TYPE(
	GeneParameters,
	/**
        \brief number of loci per gene
    */
     ( (Int64)(locusNumber) )
    
    /**
        \brief number of alleles per locus
     */
     ( (Array<Double>)(alleleNumber) )
    
	/**
		\brief Transmissibility of genes.
		
		If only one entry present, used for all genes.
	*/
	( (Array<Double>)(transmissibility) )
	
	/**
		\brief Immunity loss rate of genes.
		
		If only one entry present, used for all genes.
	*/
	( (Array<Double>)(immunityLossRate) )
	
	/**
		\brief Clinical immunity loss rate of genes.
		
		If only one entry present, used for all genes.
	*/
	( (Array<Double>)(clinicalImmunityLossRate) )
	
	/**
		\brief Relative probabilities of transitions between different genes.
		
		`mutationWeights[i][j]` is the relative probability that gene `i` will
		transition to gene `j`, given a mutation event.
		
		`mutationWeights` is normalized so that
		
		`sum(mutationWeights[i][...])` = 1.

	( (Array<Array<Double>>)(mutationWeights) )
    */
    /**
      \brief Relative probabilities of mutation rates of each locus within a gene
     */
     ( (Array<Double>)(mutationWeights) )

    /**
     \brief whether to include microsatellites
     */
     ( (Bool) (includeMicrosat))

    /**
    \brief how many microsat
    */
    ( (Int64) (microsatNumber))

    /**
    \brief number of alleles per microsat, can be different per microsat
    */
    ( (Array<Double>) (microsatAlleles))
)

/**
 \brief Type defining parameters governing intervention.
 */
ZPPJSON_DEFINE_TYPE(
                    InterventionParameters,
    /**
     \brief time when intervention starts
     */
     ((Double) (TimeStart))
                    
     /**
     \brief duration of intervention
     */
     ((Double) (duration))
    
     /**
     \brief amplitude change of biting rate. If amplitude is a, biting rate if b, then the new biting rate during intervention is b*a
     */
     ((Double) (amplitude))
     /**
     \brief magtitude of change in migration rate, if so migration is reduced to current rate * magnitude
     */
     ((Double)(IRSMRateAmplitude))
     /**
     \brief whether turn on the intervention mode
     */
     ((Bool) (includeIntervention))

)

/**
 \brief Type defining parameters governing intervention.
 */
ZPPJSON_DEFINE_TYPE(
                    MDAParameters,
                    /**
                     \brief time when MDA starts
                     */
                    ((Double) (TimeStartMDA))
                    
                    /**
                     \brief time interval between two MDA
                     */
                    ((Double) (interval))

                    /**
                     \brief total number times implementing MDA
                     */
                    ((Int64) (totalNumber))

                    /**
                     \brief failure rate per strain
                     */
                    ((Double) (strainFailRate))
                    
                    /**
                     \brief failure rate per host
                     */
                    ((Double) (hostFailRate))

                    /**
                     \brief drug effective length
                     */
                    ((Double) (drugEffDuration))
                    /**
                     \brief amplitude of change in migration rate, if so migration is reduced to current rate * amplitude
                     */
                     ((Double)(MDAMRateChange))
                    /**
                     \brief whether turn on the MDA mode
                     */
                    ((Bool) (includeMDA))
                    
                    )

ZPPJSON_DEFINE_TYPE(
     HostFollowingParameters,
     /**
                     \brief number of hosts to sample
                     */
                    ((Int64) (HostNumber))
                    
                    /**
                     \brief whether turn on the host Sampling mode
                     */
                    ((Bool) (includeHostFollowing))
                    
                    /**
                     \brief how long to follow each host
                     */
                    ((Double) (followDuration))
                    
                    )

/**
	\brief All simulation parameters.
*/
ZPPJSON_DEFINE_TYPE(
	SimParameters,
	
	/**
		\brief Output filename of database.
		
		Absolute pathname or taken relative to the working directory the
		program was executed in.
	*/
	( (String)(dbFilename) )
	
	/**
		\brief Whether or not to overwrite database if present
	*/
	( (Bool)(overwriteDatabase) )
	
	/**
		\brief How often, in simulation time units, to commit the database.
		
		Too-frequent commits can cause database maintenance to become a
		bottleneck. If the simulation is running oddly slow, try increasing
		this number.
	*/
	( (Double)(dbCommitPeriod) )
	
	/**
		\brief Random seed for simulation.
		
		If missing or set to zero, a seed will be generated and will appear
		in the parameter values inserted into
	*/
	( (Int64)(randomSeed) )
	
    /**
     \brief selection mode: 1.specific immunity, 2. general immunity, 3. neutral
     */
     ((Int64)(selectionMode))
    
	/**
		\brief Simulation end time
	*/
	( (Double)(tEnd) )

    /**
       \brief Burnin time, all sampling starts after this burnin
    */
    ( (Double)(burnIn) )

	/**
		\brief How often to update seasonal rates.
	*/
	( (Double)(seasonalUpdateEvery) )
	
	/**
		\brief Whether or not to write out all hosts to database
	*/
	( (Bool)(outputHosts) )
	
	/**
		\brief Whether or not to write out all genes to database
	*/
	( (Bool)(outputGenes) )
    /**
     \brief Whether or not to write out all loci profiles to database
    */
    ( (Bool)(outputLoci) )
	
	/**
		\brief Whether or not to write out all strains to database
	*/
	( (Bool)(outputStrains) )
	
	/**
		\brief How often to sample hosts
	*/
	( (Double)(sampleHostsEvery) )
	
	/**
		\brief How often to sample a transmission event, in number of transmission events.
	*/
	( (Int64)(sampleTransmissionEventEvery) )
	
	/**
		\brief Number of genes in the global pool.
	*/
	( (Int64)(genePoolSize) )
	
	/**
		\brief Number of genes per pathogen strain.
	*/
	( (Int64)(genesPerStrain) )
	
	/**
		\brief Probability per gene of a mutation when a strain is picked up.
	*/
	( (Double)(pMutation) )
	
    /**
     \brief Probability mitotic gene conversion/recombination, per transmission events
    */
    ( (Double)(pIntraRecomb) )
    
    /**
     \brief percentage of gene conversion in mitotic recombination
    */
    ( (Double)(percConversion))
                   
	/**
		\brief Probability that a transmitted strain will be a recombinant.
	*/
	( (Double)(pRecombinant) )
	
	/**
		\brief Duration of liver stage (pre-expression)
	*/
	( (Double)(tLiverStage) )
	
	/**
		\brief Parameter controlling distance function.
	*/
	( (Double)(distancePower) )
	
	/**
		\brief Host-lifetime distribution, specified as a discrete PDF,
		with uniform density within each discrete chunk
	*/
	( (DiscretizedPDF)(hostLifetimeDistribution) )
	
	/**
		\brief Array of population parameters, one set for each population.
	*/
	( (Array<PopulationParameters>)(populations) )
	
	/**
		\brief Whether or not coinfection reduces transmission.
	*/
	( (Bool)(coinfectionReducesTransmission) )
	
	/**
		\brief Parameters governing within-host dynamics
		(see WithinHostParameters class).
	*/
	( (WithinHostParameters)(withinHost) )
	
	/**
		\brief Parameters governing genes (see GeneParameters class).
	*/
	( (GeneParameters)(genes) )
	
    /**
    \brief Parameters governing interventions (see InterventionParameters class).
    */
    ( (InterventionParameters)(intervention) )
    
    /**
    \brief Parameters governing MDA measures (see MDAParameters class).
    */
    ( (MDAParameters)(MDA) )
      
    /**
     \brief Parameters govering specific host following Parameters (see HostFollowingParameters class).
     */
     ((HostFollowingParameters)(following))
                    
	/**
		\brief Whether or not clinical immunity is tracked.
	*/
	( (Bool)(trackClinicalImmunity) )
                    
    /**
     \brief Probability of microsatellites mutation, per gene, per day
     */
     ((Double)(pMsMutate))
)

#endif // #ifndef DOXYGEN_SHOULD_SKIP_THIS

#endif /* defined(__malariamodel__SimParameters__) */
