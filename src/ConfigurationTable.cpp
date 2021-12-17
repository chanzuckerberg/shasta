
#include "ConfigurationTable.hpp"

namespace shasta {
   const std::vector< pair<string, string> > configurationTable = {
    {"Nanopore-Dec2019", R"zzz(# This file contains Shasta options that, as of December 2019,
# are known to work with Oxford Nanopore reads under the following 
# circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.0.5 base caller. Also known to work with
#   reads from other Guppy releases 3.0.x and 3.1.x.

# To use this configuration file, specify Shasta option "--config PathToThisFile". 
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options require root access.

# For detailed information on all available options see here:
# https://chanzuckerberg.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://chanzuckerberg.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://chanzuckerberg.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://chanzuckerberg.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 10000

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
minAlignedFraction = 0.4

[Assembly]
consensusCaller = Bayesian:guppy-3.0.5-a

)zzz"},
    {"Nanopore-UL-Dec2019", R"zzz(# This file contains Shasta options that, as of December 2019,
# are known to work with Ultra-Long (UL) Oxford Nanopore reads 
# under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.0.5 base caller. Also known to work with
#   reads from other Guppy releases 3.0.x and 3.1.x.

# To use this configuration file, specify Shasta option "--config PathToThisFile". 
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options require root access.

# For detailed information on all available options see here:
# https://chanzuckerberg.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://chanzuckerberg.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://chanzuckerberg.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://chanzuckerberg.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 50000

[MinHash]
minBucketSize = 5
maxBucketSize = 40
minFrequency = 10

[Align]
maxSkip = 60
maxDrift = 60
minAlignedMarkerCount = 400

[Assembly]
consensusCaller = Bayesian:guppy-3.0.5-a

)zzz"},
    {"Nanopore-Jun2020", R"zzz(# DO NOT USE THIS FILE IF YOU HAVE READS CREATED BY A
# GUPPY VERSION OLDER THAN 3.6.0.

# This file contains Shasta options that, as of June 2020,
# are known to work with Oxford Nanopore reads under the following 
# circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.6.0 base caller. If you have reads
#   from an older version of Guppy, use configuration file
#   Nanopore-Dec2019.conf instead.

# To use this configuration file, specify Shasta option 
# "--config AbsolutePathToThisFile". 
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options requires root access.

# For detailed information on all available options see here:
# https://chanzuckerberg.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://chanzuckerberg.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://chanzuckerberg.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://chanzuckerberg.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 10000
noCache = True

[Kmers]
# Due to the higher accuracy of Guppy 3.6.0 we use longer
# markers than usual.
k = 14

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
minAlignedFraction = 0.55
minAlignedMarkerCount = 400
sameChannelReadAlignment.suppressDeltaThreshold = 30

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
refineThreshold = 6
crossEdgeCoverageThreshold = 3

[Assembly]
consensusCaller = Bayesian:guppy-3.6.0-a
detangleMethod = 1


)zzz"},
    {"Nanopore-UL-Jun2020", R"zzz(# DO NOT USE THIS FILE IF YOU HAVE READS CREATED BY A
# GUPPY VERSION OLDER THAN 3.6.0.

# This file contains Shasta options that, as of June 2020,
# are known to work with Ultra-Long (UL) Oxford Nanopore reads 
# under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.6.0 base caller. If you have reads
#   from an older version of Guppy, use configuration file
#   Nanopore-UL-Dec2019.conf instead.

# To use this configuration file, specify Shasta option 
# "--config AbsolutePathToThisFile". 
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options requires root access.

# For detailed information on all available options see here:
# https://chanzuckerberg.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://chanzuckerberg.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://chanzuckerberg.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://chanzuckerberg.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 50000
noCache = True

[Kmers]
# Due to the higher accuracy of Guppy 3.6.0 we use longer
# markers than usual.
k = 14

[MinHash]
minBucketSize = 10
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
minAlignedFraction = 0.55
minAlignedMarkerCount = 600
sameChannelReadAlignment.suppressDeltaThreshold = 30

[ReadGraph]
maxAlignmentCount = 12


[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
refineThreshold = 6
crossEdgeCoverageThreshold = 3

[Assembly]
consensusCaller = Bayesian:guppy-3.6.0-a
detangleMethod = 1


)zzz"},
    {"Nanopore-Sep2020", R"zzz(# DO NOT USE THIS FILE IF YOU HAVE READS CREATED BY A
# GUPPY VERSION OLDER THAN 3.6.0.

# This file contains Shasta options which attempt to partially automate
# parameter selection. It is based on an earlier config, which, as of Jun 2020,
# was known to work with Oxford Nanopore reads under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.6.0 base caller. If you have reads
#   from an older version of Guppy, use configuration file
#   Nanopore-Dec2019.conf instead.

# The automation provided by this config is particularly applicable to
# low coverage or non-human samples. It also matches or exceeds continuity
# in human samples, relative to the appropriately chosen 3.6.0 or 3.6.0-UL conf.
# Automation can also be activated with parameters designed for earlier basecallers,
# if needed, but updating to guppy 3.6.0 or higher will greatly improve assembly
# quality and is therefore strongly recommended.

# To use this configuration file, specify Shasta option 
# "--config AbsolutePathToThisFile". 
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options requires root access.

# For detailed information on all available options see here:
# https://chanzuckerberg.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://chanzuckerberg.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://chanzuckerberg.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://chanzuckerberg.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 10000
noCache = True

[Kmers]
# Due to the higher accuracy of Guppy 3.6.0 we use longer
# markers than usual.
k = 14

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This method uses the observed distribution of alignment stats to choose a cutoff for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Automatically determine this using PeakFinder
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-3.6.0-a
detangleMethod = 2


)zzz"},
    {"Nanopore-UL-Sep2020", R"zzz(# DO NOT USE THIS FILE IF YOU HAVE READS CREATED BY A
# GUPPY VERSION OLDER THAN 3.6.0.

# This file contains Shasta options which attempt to partially automate
# parameter selection. It is based on an earlier config, which, as of Jun 2020,
# was known to work with Oxford Nanopore reads under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.6.0 base caller. If you have reads
#   from an older version of Guppy, use configuration file
#   Nanopore-UL-Dec2019.conf instead.

# The automation provided by this config is particularly applicable to
# low coverage or non-human samples. It also matches or exceeds continuity
# in human samples, relative to the appropriately chosen 3.6.0 or 3.6.0-UL conf.
# Automation can be activated with parameters designed for earlier basecallers,
# if needed, but updating to guppy 3.6.0 or higher will greatly improve assembly
# quality and is therefore strongly recommended.

# In most cases, for best performance on a large assembly 
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options requires root access.

# For detailed information on all available options see here:
# https://chanzuckerberg.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://chanzuckerberg.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://chanzuckerberg.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://chanzuckerberg.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 50000
noCache = True

[Kmers]
# Due to the higher accuracy of Guppy 3.6.0 we use longer
# markers than usual.
k = 14

[MinHash]
minBucketSize = 10
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This method uses the observed distribution of alignment stats to choose a cutoff for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Automatically determine this using PeakFinder
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-3.6.0-a
detangleMethod = 2


)zzz"},
    {"Nanopore-UL-iterative-Sep2020", R"zzz(# This configuration file is EXPERIMENTAL and should only
# be used under the following conditions:
# - Nanopore reads created by Guppy 3.6.0 or newer.
# - Ultra-Long (UL) reads with typical N50 80Kb or more.
# - High coverage 80X.
# Iterative assembly results in some separation of haplotypes
# and does better at resolving segmental duplications.

[Reads]
minReadLength = 30000 
noCache = True

[Kmers]
k = 10 

[MinHash]
minBucketSize = 10 
maxBucketSize = 40 
minFrequency = 5 

[Align]
alignMethod = 3 
matchScore = 6
gapScore = -3 
downsamplingFactor = 0.05 
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1
sameChannelReadAlignment.suppressDeltaThreshold = 30 

[ReadGraph]
maxAlignmentCount = 12
creationMethod = 2

[MarkerGraph]
minCoveragePerStrand = 3
simplifyMaxLength = 10,100
crossEdgeCoverageThreshold = 3 

[Assembly]
detangleMethod = 2 
consensusCaller = Bayesian:guppy-3.6.0-a
iterative = True

)zzz"},
    {"Nanopore-OldGuppy-Sep2020", R"zzz(# This file contains Shasta options which attempt to partially automate
# parameter selection. It is based on an earlier config, which, as of Jun 2020,
# was known to work with Oxford Nanopore reads under the following circumstances:

# - Human genome assembly.
# - Coverage between 40x and 80x. If you have more coverage than that,
#   you can use option "--Reads.minReadLength" to adjust coverage as desired.
# - Reads from Guppy 3.0.5 base caller. Also known to work with
#   reads from other Guppy releases 3.0.x and 3.1.x.

# The automation provided by this config is particularly applicable to
# low coverage or non-human samples. It also matches or exceeds continuity
# in human samples, relative to the appropriately chosen config file.
# Updating to guppy 3.6.0 or higher will greatly improve assembly
# quality and is therefore strongly recommended.

# To use this configuration file, specify Shasta option
# "--config AbsolutePathToThisFile".
# If you specify any conflicting values on the command line,
# the values specified on the command line take precedence.

# In most cases, for best performance on a large assembly
# you will usually also want to use the following options, which 
# cannot be specified in a configuration file:
# --memoryMode filesystem
# --memoryBacking 2M
# Using these options requires root access.

# For detailed information on all available options see here:
# https://chanzuckerberg.github.io/shasta/CommandLineOptions.html

# For information on running a small assembly for which 
# performance is not essential see here:
# https://chanzuckerberg.github.io/shasta/QuickStart.html

# For more information on running an assembly see here:
# https://chanzuckerberg.github.io/shasta/Running.html

# For information on optimizing assembly performance see here:
# https://chanzuckerberg.github.io/shasta/Performance.html



[Reads]
# If you have extra coverage, use this option to adjust coverage.
minReadLength = 10000
noCache = True

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This method uses the observed distribution of alignment stats to choose a cutoff for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Automatically determine this using PeakFinder
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-3.0.5-a
detangleMethod = 2

)zzz"},
    {"Nanopore-Plants-Apr2021", R"zzz(# This was used in April 2021 for an assembly of Lolium Perenne
# under the following conditions:

# - Oxford Nanopore reads, basecalled with Guppy 4.0.14.
# - Coverage: 30x.
# - Reads N50: 62 Kb. 

# The assembly was run:

# - Using Shasta at commit 266d4c8e65ff408db2dfd1381bb14648be83dad3.
# - On AWS Graviton2 (ARM) r6g.16xlarge instance type, 512 GB, 64 vCPUs.
# - Using memory options --memoryMode filesystem --memoryBacking 2M.

# Assembly results:

# - Assembled 2.185 Mb of sequence. Estimated genome size 
#   is 2.46 Gb from k-mer analysis, 2.72 Gb from flow cytometry.
# - Assembly N50: 5.5 Mb.
# - Elapsed time for assembly: 95 minutes.

# Many thanks to Dario Copetti (Molecular Plant Breeding, ETH Zurich, Switzerland)
# for providing access to the reads in advance of publication. 

# Also see Shasta issue #200 for some discussion
# https://github.com/chanzuckerberg/shasta/issues/200



[Reads]
noCache = True

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minHashIterationCount = 50
minFrequency = 5

[Align]
downsamplingFactor = 0.05
sameChannelReadAlignment.suppressDeltaThreshold = 30
maxSkip = 60
maxDrift = 20
maxTrim = 60
minAlignedMarkerCount = 200
minAlignedFraction = 0.3

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-3.6.0-a
detangleMethod = 2


)zzz"},
    {"Nanopore-Oct2021", R"zzz(# This is known to work at least under the following conditions:
# - Oxford Nanopore reads.
# - Guppy 5 base caller.
# - Human genome.
# - Coverage 30x to 80x.

[Reads]
minReadLength = 10000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
# (ReadGraph.creationMethod 2).
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This uses the observed distribution of alignment statistics to choose thresholds for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-a
detangleMethod = 2


)zzz"},
    {"Nanopore-UL-Oct2021", R"zzz(# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long (UL) reads.
# - Guppy 5 base caller.
# - Human genome.
# - Coverage 30x to 80x.



[Reads]
minReadLength = 50000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 10
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
# (ReadGraph.creationMethod 2).
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This uses the observed distribution of alignment statistics to choose thresholds for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-a
detangleMethod = 2


)zzz"},
    {"HiFi-Oct2021", R"zzz(# This was tested on the  "HG002 Data Freeze (v1.0) 
# Recommended downsampled data mix", a 34x HiFi dataset
# for HG002. See here for more information:
# https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0#hg002-data-freeze-v10-recommended-downsampled-data-mix
# It assembled 3078 Mb of sequence with an N50 of 35 Mb.
# The quality of assembled sequence was estimated
# at Q = 46 dB via k-mer based methods. It is limited
# by the fact that this creates a haploid assembly
# in which one of the alleles of each heterozygous locus
# is removed.



[Reads]
minReadLength = 8000
noCache = True

[Kmers]
k = 14

[MinHash]
hashFraction = 0.05
minHashIterationCount = 100
minFrequency = 3
minBucketSize = 10
maxBucketSize = 60

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
minAlignedFraction = 0.97
minAlignedMarkerCount = 200
maxSkip = 6
maxDrift = 4
maxTrim = 2

[ReadGraph]
maxAlignmentCount = 30
maxChimericReadDistance = 2

[MarkerGraph]
minCoverage = 6
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

[Assembly]
consensusCaller = Modal
detangleMethod = 2


)zzz"}
    };
}
